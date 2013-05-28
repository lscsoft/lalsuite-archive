#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
Text-mode progress bars
"""
__copyright__ = "Copyright 2010, Leo Singer"
__author__    = "Leo Singer <leo.singer@ligo.org>"
__all__       = ["ProgressBar"]


# TODO: Implement the context manager API (but wait until clusters have Python >= 2.5).
# TODO: Is it better to print the progress bar to stderr or stdout?


import os
import sys


# From http://stackoverflow.com/questions/566746/how-to-get-console-window-width-in-python
def getTerminalSize():
    """
    returns (lines:int, cols:int)
    """
    import os, struct
    def ioctl_GWINSZ(fd):
        import fcntl, termios
        return struct.unpack("hh", fcntl.ioctl(fd, termios.TIOCGWINSZ, "1234"))
    # try stdin, stdout, stderr
    for fd in (0, 1, 2):
        try:
            return ioctl_GWINSZ(fd)
        except:
            pass
    # try os.ctermid()
    try:
        fd = os.open(os.ctermid(), os.O_RDONLY)
        try:
            return ioctl_GWINSZ(fd)
        finally:
            os.close(fd)
    except:
        pass
    # try `stty size`
    try:
        return tuple(int(x) for x in os.popen("stty size", "r").read().split())
    except:
        pass
    # try environment variables
    try:
        return tuple(int(os.getenv(var)) for var in ("LINES", "COLUMNS"))
    except:
        pass
    # i give up. return default.
    return (25, 80)


class ProgressBar:
    """Display a text progress bar."""

    def __init__(self, text='Working', max=1, value=0, textwidth=24, fid=None):
        if fid is None:
            self.fid = sys.stderr
        self.isatty = os.isatty(self.fid.fileno())
        self.text = text
        self.max = max
        self.value = value
        self.textwidth = textwidth
        self.twiddle = 0

    def iterate(self, iterable, format="%s"):
        """Use as a target of a for-loop to issue a progress update for every
        iteration. For example:

        progress = ProgressBar()
        for text in progress.iterate(["foo", "bar", "bat"]):
            ...
        """

        # If iterable has a definite length, then set the maximum value of the
        # progress bar. Else, set the maximum value to -1 so that the progress
        # bar displays indeterminate progress (scrolling dots).
        try:
            length = len(iterable)
        except TypeError:
            self.max = -1
        else:
            self.max = length

        # Iterate over the input, updating the progress bar for each element.
        for i, item in enumerate(iterable):
            self.update(i, format % item)
            yield item

    def show(self):
        """Redraw the text progress bar."""
        from math import floor, ceil

        if len(self.text) > self.textwidth:
            label = self.text[0:self.textwidth]
        else:
            label = self.text.rjust(self.textwidth)

        terminalSize = getTerminalSize()
        if terminalSize is None:
            terminalSize = 80
        else:
            terminalSize = terminalSize[1]

        barWidth = terminalSize - self.textwidth - 10

        if self.value is None or self.value < 0:
            if self.twiddle == 0:
                pattern = ' ..'
                self.twiddle = 1
            elif self.twiddle == 1:
                pattern = '. .'
                self.twiddle = 2
            else:
                pattern = '.. '
                self.twiddle = 0
            barSymbols = (pattern * int(ceil(barWidth/3.0)))[0:barWidth]
            progressFractionText = '   . %'
        else:
            solidBlock = '|'
            partialBlocks = ' .:!'
            blank = ' '

            progressFraction = float(self.value) / self.max
            nBlocks = progressFraction * barWidth
            nBlocksInt = int(floor(progressFraction * barWidth))
            partialBlock = partialBlocks[int(floor((nBlocks - nBlocksInt) * len(partialBlocks)))]
            nBlanks = barWidth - nBlocksInt - 1
            barSymbols = (solidBlock * nBlocksInt) + partialBlock + (blank * nBlanks)
            progressFractionText = ('%.1f%%' % (100*progressFraction)).rjust(6)

        print >>self.fid, '\r\x1B[1m' + label + '\x1B[0m [' + barSymbols + ']' + progressFractionText,
        self.fid.flush()

    def update(self, value = None, text = None):
        """Redraw the progress bar, optionally changing the value and text."""
        old_text = self.text
        if text is not None:
            self.text = text
        if value is not None:
            self.value = value
        if self.isatty:
            self.show()
        elif self.text != old_text:
            print >>self.fid, self.text


def demo():
    """Demonstrate progress bar."""
    from time import sleep
    maxProgress = 1000
    progressbar = ProgressBar(max=maxProgress)
    for i in range(-100,maxProgress):
        sleep(0.01)
        progressbar.update(i)
    for s in progressbar.iterate(range(maxProgress)):
        sleep(0.01)
    for s in progressbar.iterate(range(maxProgress), format='iteration %d'):
        sleep(0.01)
    print >>sys.stderr, ''


if __name__ == '__main__':
    demo()
