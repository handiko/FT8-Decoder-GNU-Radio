#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: FT8 Multi-Signal Receiver (sound card)
# Description: Multi-signal FT8 decoder from a sound card
# GNU Radio version: 3.10.10.0

from PyQt5 import Qt
from gnuradio import qtgui
from gnuradio import audio
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio import gr
from gnuradio.fft import window
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
import ft8_rx_ft8_clock as ft8_clock  # embedded python block
import ft8_rx_ft8_demod as ft8_demod  # embedded python block
import ft8_rx_ft8_ldpc as ft8_ldpc  # embedded python block
import ft8_rx_ft8_report as ft8_report  # embedded python block
import ft8_rx_ft8_search as ft8_search  # embedded python block
import ft8_rx_ft8_subtract as ft8_subtract  # embedded python block
import sip


def snipfcn_snippet_table(self):
    # Put the FT8 decode table to the right of the spectrum/waterfall
    self.top_grid_layout.addWidget(self.ft8_report.make_widget(), 0, 1, 2, 1)
    self.top_grid_layout.setColumnStretch(0, 3)
    self.top_grid_layout.setColumnStretch(1, 2)
    self.resize(1400, 800)


def snippets_main_after_init(tb):
    snipfcn_snippet_table(tb)

class ft8_rx(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "FT8 Multi-Signal Receiver (sound card)", catch_exceptions=True)
        Qt.QWidget.__init__(self)
        self.setWindowTitle("FT8 Multi-Signal Receiver (sound card)")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except BaseException as exc:
            print(f"Qt GUI: Could not set Icon: {str(exc)}", file=sys.stderr)
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "ft8_rx")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Variables
        ##################################################
        self.wav_file = wav_file = "C:/Users/me/AppData/Local/WSJT-X/save/000000_000000.wav"
        self.samp_rate = samp_rate = 48000
        self.pc_clock = pc_clock = True
        self.passes = passes = 3
        self.ft8_rate = ft8_rate = 12000
        self.f_min = f_min = 200
        self.f_max = f_max = 3000

        ##################################################
        # Blocks
        ##################################################

        self.xlate_display = filter.freq_xlating_fir_filter_fcf(3, firdes.low_pass(1.0, ft8_rate, 1900, 200), 2000, ft8_rate)
        self.waterfall = qtgui.waterfall_sink_c(
            1024, #size
            window.WIN_BLACKMAN_hARRIS, #wintype
            2000, #fc
            (ft8_rate / 3), #bw
            "Waterfall (audio Hz)", #name
            1, #number of inputs
            None # parent
        )
        self.waterfall.set_update_time(0.10)
        self.waterfall.enable_grid(False)
        self.waterfall.enable_axis_labels(True)

        self.waterfall.disable_legend()


        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        colors = [0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in range(1):
            if len(labels[i]) == 0:
                self.waterfall.set_line_label(i, "Data {0}".format(i))
            else:
                self.waterfall.set_line_label(i, labels[i])
            self.waterfall.set_color_map(i, colors[i])
            self.waterfall.set_line_alpha(i, alphas[i])

        self.waterfall.set_intensity_range(-100, -40)

        self._waterfall_win = sip.wrapinstance(self.waterfall.qwidget(), Qt.QWidget)

        self.top_grid_layout.addWidget(self._waterfall_win, 1, 0, 1, 1)
        for r in range(1, 2):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.spectrum = qtgui.freq_sink_c(
            2048, #size
            window.WIN_BLACKMAN_hARRIS, #wintype
            2000, #fc
            (ft8_rate / 3), #bw
            "Spectrum (audio Hz)", #name
            1,
            None # parent
        )
        self.spectrum.set_update_time(0.20)
        self.spectrum.set_y_axis((-120), (-20))
        self.spectrum.set_y_label('Relative Gain', 'dB')
        self.spectrum.set_trigger_mode(qtgui.TRIG_MODE_FREE, 0.0, 0, "")
        self.spectrum.enable_autoscale(False)
        self.spectrum.enable_grid(True)
        self.spectrum.set_fft_average(0.2)
        self.spectrum.enable_axis_labels(True)
        self.spectrum.enable_control_panel(False)
        self.spectrum.set_fft_window_normalized(False)

        self.spectrum.disable_legend()


        labels = ['', '', '', '', '',
            '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
            "magenta", "yellow", "dark red", "dark green", "dark blue"]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]

        for i in range(1):
            if len(labels[i]) == 0:
                self.spectrum.set_line_label(i, "Data {0}".format(i))
            else:
                self.spectrum.set_line_label(i, labels[i])
            self.spectrum.set_line_width(i, widths[i])
            self.spectrum.set_line_color(i, colors[i])
            self.spectrum.set_line_alpha(i, alphas[i])

        self._spectrum_win = sip.wrapinstance(self.spectrum.qwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._spectrum_win, 0, 0, 1, 1)
        for r in range(0, 1):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.lpf_dec4 = filter.fir_filter_fff(
            (int(samp_rate // ft8_rate)),
            firdes.low_pass(
                1,
                samp_rate,
                4500,
                1000,
                window.WIN_HAMMING,
                6.76))
        self.ft8_subtract = ft8_subtract.blk()
        self.ft8_search = ft8_search.blk(f_min=f_min, f_max=f_max, syncmin=1.2, maxcand=400)
        self.ft8_report = ft8_report.blk(passes=passes, console=True)
        self.ft8_ldpc = ft8_ldpc.blk(max_iter=40)
        self.ft8_demod = ft8_demod.blk()
        self.ft8_clock = ft8_clock.blk(samp_rate=ft8_rate, decode_at=14.4, pc_clock=pc_clock, clock_offset=0.0)
        self.audio_source_0 = audio.source(samp_rate, '', True)


        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.ft8_clock, 'slot'), (self.ft8_search, 'slot'))
        self.msg_connect((self.ft8_demod, 'llr'), (self.ft8_ldpc, 'llr'))
        self.msg_connect((self.ft8_demod, 'audio'), (self.ft8_subtract, 'audio'))
        self.msg_connect((self.ft8_ldpc, 'cw'), (self.ft8_report, 'cw'))
        self.msg_connect((self.ft8_report, 'new'), (self.ft8_subtract, 'new'))
        self.msg_connect((self.ft8_search, 'cand'), (self.ft8_demod, 'cand'))
        self.msg_connect((self.ft8_subtract, 'slot'), (self.ft8_search, 'slot'))
        self.connect((self.audio_source_0, 0), (self.lpf_dec4, 0))
        self.connect((self.lpf_dec4, 0), (self.ft8_clock, 0))
        self.connect((self.lpf_dec4, 0), (self.xlate_display, 0))
        self.connect((self.xlate_display, 0), (self.spectrum, 0))
        self.connect((self.xlate_display, 0), (self.waterfall, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "ft8_rx")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_wav_file(self):
        return self.wav_file

    def set_wav_file(self, wav_file):
        self.wav_file = wav_file

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.lpf_dec4.set_taps(firdes.low_pass(1, self.samp_rate, 4500, 1000, window.WIN_HAMMING, 6.76))

    def get_pc_clock(self):
        return self.pc_clock

    def set_pc_clock(self, pc_clock):
        self.pc_clock = pc_clock
        self.ft8_clock.pc_clock = self.pc_clock

    def get_passes(self):
        return self.passes

    def set_passes(self, passes):
        self.passes = passes
        self.ft8_report.passes = self.passes

    def get_ft8_rate(self):
        return self.ft8_rate

    def set_ft8_rate(self, ft8_rate):
        self.ft8_rate = ft8_rate
        self.xlate_display.set_taps(firdes.low_pass(1.0, self.ft8_rate, 1900, 200))
        self.spectrum.set_frequency_range(2000, (self.ft8_rate / 3))
        self.waterfall.set_frequency_range(2000, (self.ft8_rate / 3))

    def get_f_min(self):
        return self.f_min

    def set_f_min(self, f_min):
        self.f_min = f_min
        self.ft8_search.f_min = self.f_min

    def get_f_max(self):
        return self.f_max

    def set_f_max(self, f_max):
        self.f_max = f_max
        self.ft8_search.f_max = self.f_max




def main(top_block_cls=ft8_rx, options=None):

    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    snippets_main_after_init(tb)
    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    qapp.exec_()

if __name__ == '__main__':
    main()
