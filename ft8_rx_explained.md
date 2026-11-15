# FT8 Multi-Signal Receiver (`ft8_rx.grc`, `ft8_rx_staged.grc`)

There are two flowgraphs. Both have the same front end and display, decode identically, and run at almost the same speed. `ft8_rx.grc` puts the whole decoder in one block. `ft8_rx_staged.grc` splits it into six visible stages (see "Staged version" below). The plain-English story of the staged version is `ft8_rx_plain_english.md`.

Sound-card FT8 decoder for GNU Radio 3.10 on Windows (radioconda). It uses stock blocks plus one Embedded Python Block (numpy only, no OOT modules, no DLLs). Every 15 s it prints all decoded messages as a WSJT-X style table (`UTC dB DT Freq Message`), both in the GRC console and in the QT window, next to a spectrum and a waterfall.

## Flowgraph

```
Audio Source (48 kHz) → Low Pass Filter (4.5 kHz, ÷4) → 12 kHz ─┬→ FT8 Decoder (EPB) → console + table
                                                               └→ Freq Xlating FIR (2 kHz, ÷3) → Spectrum + Waterfall (0–4 kHz)
Python Snippet: puts the decoder's table into the QT window
Disabled test path: Wav File Source → Throttle (for WSJT-X .wav files)
```

| Operation | Block | Why |
|---|---|---|
| Sound-card capture | Audio Source | stock |
| Anti-alias + decimate to 12 kS/s | Low Pass Filter | stock |
| Display band 0–4 kHz | Freq Xlating FIR → QT Frequency / Waterfall Sink | stock; the shift centres the plots on 2 kHz so the axis reads 0–4000 Hz audio |
| Slot timing, Costas search, demod, LDPC, CRC, unpack, subtraction | EPB `FT8 Decoder` | no stock blocks for these; gr-fec's LDPC was rejected because the subtract-and-decode-again loop needs the decode result inside the same block |
| Table in the QT window | Python Snippet (Main – After Init) | GRC has no stock text table widget |

## Decoder (inside the EPB)

1. **Slot timing.** Slots follow the PC clock (:00/:15/:30/:45). The block ties sample counts to the PC clock and re-anchors every slot, which absorbs sound-card clock drift (tested with an 80 ppm error and random delivery jitter). Decoding starts 14.4 s into the slot in a worker thread, so the audio stream never stalls. A slot joined with less than 13 s of audio is skipped.
2. **Coarse search.** Power spectrogram (3.125 Hz × 40 ms) and a 7×7 Costas score for every frequency and lag from −2.5 s to +2.5 s; up to 400 candidates.
3. **Per candidate.** One FFT of the whole slot, then a small IFFT per candidate: shift to 0 Hz and decimate to 200 S/s. Fine time (±50 ms) and frequency (±2.5 Hz) sync, then a 32-point FFT per symbol → 8 tone amplitudes.
4. **Soft bits.** 1-, 2- and 3-symbol coherent metrics plus a per-bit normalized metric, the same four variants WSJT-X tries.
5. **LDPC(174,91).** Sum-product BP with all candidates × 4 metrics in one numpy batch, then the CRC-14 check, at most 36 hard errors, and a valid unpack.
6. **Subtraction (passes = 3).** Every decoded signal is rebuilt (GFSK, BT = 2), fitted for amplitude and phase, and subtracted. The residual is searched again, which finds weak signals hidden under strong ones.
7. **Messages.** Standard, /R and /P, nonstandard calls with hashes (`<CALL>` or `<...>`), free text, DXpedition, Field Day, RTTY Roundup, EU VHF and telemetry. The callsign hash table persists across slots.

## Measured on 22 real recorded slots (ft8_lib test set: 9 saved by WSJT-X, 13 from WebSDR receivers)

| Decoder | Messages decoded |
|---|---|
| WSJT-X (reference list) | 353 |
| **This decoder, passes = 3** | **316 of those (90 %)**, plus 16 valid messages missing from the WSJT-X list (real calls with matching partners; two also decoded by ft8_lib) |
| This decoder, passes = 1 | 268 (76 %) |
| ft8_lib decode_ft8 | 259 (73 %) |

- False decodes in 30 slots of pure noise: 0.
- dB column vs WSJT-X: mean offset 0 dB, spread ±3 dB.
- DT matches WSJT-X to 0.1 s.
- Decode time: about 1.8 s per slot on average with passes = 3 (0.6 s with passes = 1), up to about 2.5 s for the busiest slots, on the test machine (one core, numpy).

## Use

**Live:**
1. Feed the radio's line/data/speaker output to Line In.
2. Set USB mode with a 3 kHz filter and the radio's NR, notch and narrow filters off. Turn off Windows "audio enhancements" and mic boost, and keep the level well below clipping.
3. Leave the Audio Source device name blank to use the default recording device, or enter the device name.
4. Run (F6). The first table appears about 17 s after the first full slot.

**Test with WSJT-X recordings:**
1. In WSJT-X, choose Save → Save all; the recordings are 12 kHz mono .wav files.
2. Disable Audio Source and Low Pass Filter, and enable Wav File Source and Throttle.
3. Set `wav_file`, and set `pc_clock = False` so slots are counted from the start of the file.

**Variables:**
- `f_min` / `f_max`: search range in Hz (default 200–3000).
- `passes`: 1 is fastest; 3 gives the most decodes.
- The EPB parameter `clock_offset` adds a fixed offset (in seconds) to the PC clock if DT is consistently off; the time-sync software already handles the clock itself.

The `decodes` message port gives one PMT dict per decode (`utc`, `snr`, `dt`, `freq`, `msg`), so it can feed a Message Debug block or a UDP sink later.

## Staged version (`ft8_rx_staged.grc`)

This version uses the same decoder, split into six Embedded Python Blocks that pass one message per slot and per pass. Each block does one job, and the subtraction loop shows up in GRC as a message connection going backwards.

```
                    ┌──────────────────────── residual audio (passes 2, 3) ─────────────────────┐
                    ▼                                                                            │
 LPF ─▶ 1 Slot Clock ─slot─▶ 2 Sync Search ─cand─▶ 3 Demodulator ─llr─▶ 4 LDPC + CRC ─cw─▶ 5 Unpack + Report ─▶ table / console / 'decodes'
                                                         │ audio                              │ new decodes
                                                         └──────────────▶ 6 Subtract ◀────────┘
```

| Stage | Input → output | What it does |
|---|---|---|
| 1 Slot Clock | 12 kHz stream → `slot` | Cuts the audio into 15 s UTC slots using the PC clock and follows sound-card drift. File mode counts slots from the first sample. |
| 2 Sync Search | `slot` → `cand` | Spectrogram + Costas 7×7 score over frequency and lag → candidate list. Normalizes the audio on pass 0. |
| 3 Demodulator | `cand` → `llr`, `audio` | Shift to 0 Hz, decimate to 200 S/s, fine sync, 79×8 tone amplitudes, 4 soft-bit metrics. Drops candidates with ≤ 6 of 21 sync tones correct. |
| 4 LDPC + CRC | `llr` → `cw` | Batch BP decoding, CRC-14, hard-error limit. |
| 5 Unpack + Report | `cw` → `new`, `decodes` | Removes duplicates, unpacks, measures SNR/DT/Freq. Decides whether another pass runs (`passes`). Prints the table and feeds the QT table. |
| 6 Subtract | `audio` + `new` → `slot` | Rebuilds each new decode as GFSK, fits it, subtracts it, and sends the residual to stage 2 as the next pass. |

**Cost of splitting, measured on the same 22 slots:**
- Decodes: identical (316 matched + 16 extra).
- Time: 1.78 s per slot against 1.69 s for the single block, about +5% from converting messages; each 180 000-sample slot costs about 15 ms per hop.
- Running inside GNU Radio with the QT sinks: 1.6–2.0 s per slot.
- GNU Radio runs each block's message handler in its own thread, so stages of consecutive slots can overlap.

**Messages:** every message is a PMT dict carrying `slot`, `label`, `pretty`, `pass` and `t0` plus that stage's arrays (`audio`, `f0`/`xdt`, `llr`, `cw`, `cs`, `base`, …). To look inside, connect a Message Debug block to any arrow. The audio-carrying messages are large, so for those print only the field names.

**Note for file tests:** with Wav File Source `repeat = False`, GNU Radio stops the message blocks when the file ends, so the last slot may not print. Keep `repeat = True` (the default here).

## Limits

- No a-priori decoding and no OSD fallback. These are what give WSJT-X its last 1–3 dB.
- Only FT8; FT4 is not supported.
