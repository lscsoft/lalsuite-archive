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
__author__ = "Leo Singer <leo.singer@ligo.org>"
__all__ = ["ProgressBar"]


import math
import os
import struct
import sys


# From http://stackoverflow.com/questions/566746
def getTerminalSize():
    """
    returns (lines:int, cols:int)
    """

    def ioctl_GWINSZ(fd):
        # These two imports are only present on POSIX systems, so they must be
        # guarded by a try block.
        import fcntl
        import termios
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
    # try environment variables
    try:
        return tuple(int(os.getenv(var)) for var in ("LINES", "COLUMNS"))
    except:
        pass
    # i give up. return default.
    return (25, 80)


class ProgressBar:
    """Display a text progress bar.

    A final line feed is printed when the ProgressBar is garbage collected.
    Explicitly deleting the object can force a line feed when desired.  As
    an alternative, using the ProgressBar as a context manager will ensure
    a final line feed is printed when the code block within which the
    ProgressBar is being used exits.

    Example:

    >>> with ProgressBar(max=3) as pb:
    ...    pb.update(1)
    ...    pb.update(2)
    ...    pb.update(3)
    ...
    """

    def __init__(
            self, text='Working', max=1, value=0, textwidth=24, fid=None,
            sequence=' .:!|', twiddle_sequence=(' ..', '. .', '.. ')):
        if fid is None:
            self.fid = sys.stderr
        self.isatty = os.isatty(self.fid.fileno())
        self.text = text
        self.max = max
        self.value = value
        self.textwidth = textwidth
        self.sequence = sequence
        self.twiddle_sequence = twiddle_sequence
        self.twiddle = 0
        self.linefed = False

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
            pattern = self.twiddle_sequence[
                self.twiddle % len(self.twiddle_sequence)]
            self.twiddle += 1
            barSymbols = (pattern * int(math.ceil(barWidth/3.0)))[0:barWidth]
            progressFractionText = '   . %'
        else:
            progressFraction = float(self.value) / self.max

            nBlocksFrac, nBlocksInt = math.modf(
                max(0.0, min(1.0, progressFraction)) * barWidth)
            nBlocksInt = int(nBlocksInt)

            partialBlock = self.sequence[
                int(math.floor(nBlocksFrac * len(self.sequence)))]

            nBlanks = barWidth - nBlocksInt - 1
            barSymbols = (self.sequence[-1] * nBlocksInt) + partialBlock + \
                (self.sequence[0] * nBlanks)
            barSymbols = barSymbols[:barWidth]
            progressFractionText = ('%.1f%%' % (100*progressFraction)).rjust(6)

        print >>self.fid, '\r\x1B[1m' + label + '\x1B[0m [' + barSymbols + \
            ']' + progressFractionText,
        self.fid.flush()
        self.linefed = False

    def update(self, value=None, text=None):
        """Redraw the progress bar, optionally changing the value and text
        and return the (possibly new) value.  For I/O performance, the
        progress bar might not be written to the terminal if the text does
        not change and the value changes by too little.  Use .show() to
        force a redraw."""
        redraw = False
        if text is not None:
            redraw = text != self.text
            self.text = text
        if value is not None:
            redraw |= self.max == 0 or round(value/(0.0003*self.max)) != \
                round(self.value/(0.0003*self.max))
            self.value = value
        if redraw:
            if self.isatty:
                self.show()
            else:
                print >>self.fid, self.text
        return self.value

    def increment(self, delta=1, text=None):
        """Redraw the progress bar, incrementing the value by delta
        (default=1) and optionally changing the text.  Returns the
        ProgressBar's new value.  See also .update()."""
        return self.update(value=min(self.max, self.value + delta), text=text)

    def linefeed(self):
        if not self.linefed:
            print >>self.fid
            self.fid.flush()
            self.linefed = True

    def __enter__(self):
        self.show()
        return self

    def __exit__(self, exc_type, exc_value, tb):
        try:
            self.linefeed()
        except:
            pass

    def __del__(self):
        self.linefeed()


def demo():
    """Demonstrate progress bar."""
    from time import sleep
    maxProgress = 1000
    with ProgressBar(max=maxProgress) as progressbar:
        for i in range(-100, maxProgress):
            sleep(0.01)
            progressbar.update(i)
    progressbar2 = ProgressBar(max=maxProgress)
    for s in progressbar2.iterate(range(maxProgress)):
        sleep(0.01)
    for s in progressbar2.iterate(range(maxProgress), format='iteration %d'):
        sleep(0.01)


if __name__ == '__main__':
    demo()
