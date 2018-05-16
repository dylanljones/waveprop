# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones

This module enables updating the console output using the line object below

"""
from sys import stdout as out

WIDTH = 100  # Maximum line width
_current_line = None  # current line object


def line_exists():
    """Check if line exists

    Returns
    -------
    line_exists : bool
    """
    global _current_line
    return _current_line is not None


def start(text=""):
    """Start new updatable line

    Parameters
    ----------
    text : str
        title of line
    """
    global _current_line
    _current_line = StatusLine(text, WIDTH)


def update_status(status):
    """ Update status of line

    Parameters
    ----------
    status : str
        text to display
    """
    global _current_line
    if _current_line is not None:
        _current_line.update(status)


def update_progress(i, n):
    """ Update Progress of line

    Parameters
    ----------
    i : int
        current number
    n : int
        maximal number
    """
    global _current_line
    prog = progress_string(i, n)
    if _current_line is not None:
        _current_line.update(prog)


def end(text=None):
    """End the current line

    Parameters
    ----------
    text : str
        optional end text
    """
    global _current_line
    if _current_line is not None:
        _current_line.end(text)
        _current_line = None


class StatusLine:
    """ Line object

    Object to represent a updateable line with methods to manipulate the output
    """

    def __init__(self, text, width=WIDTH):
        """ Create Line

        Parameters
        ----------
        text : str
            line title
        width : int
            maximal number of charackters
        """
        self.width = width
        self.text = text
        self.start()

    def start(self):
        """Initialize updatable line"""
        string = fix_length(self.text + "...", self.width)
        out.write("\r" + string)
        out.flush()

    def update(self, status):
        """ Update line text

        Parameters
        ----------
        status : str
            text to display after line title
        """
        string = self._build_string(status)
        out.write("\r" + string)
        out.flush()

    def end(self, text):
        """Close the current line

        Parameters
        ----------
        text : str
            text to display permanently
        """
        if text is None:
            string = self._build_string("Completed!")
        else:
            string = self._build_string(text)
        out.write("\r" + string + "\n")
        out.flush()

    def _build_string(self, text):
        """ Adds the text to the line title and cuts tail if the line is to long

        Parameters
        ----------
        text : str
            text to add

        Returns
        -------
        line_txt : str
            text to display
        """
        if self.text:
            string = self.text + ": " + text
        else:
            string = text
        return fix_length(string, self.width)


def progress_string(i, n):
    """ Helper to build progress string

    Parameters
    ----------
    i : int
        current number
    n : int
        maximal number


    Returns
    -------
    progress_str : str
    """
    width = len(str(n))
    string = "({0:{width}d}/{1:d})".format(i, n, width=width)
    return string


def fix_length(msg, length):
    """ Checks if text is to long and replaces overlay by "..."

    Parameters
    ----------
    msg : str
        text to check
    length : int
        maximal number of allowed charackters

    Returns
    -------
    line : str
    """
    string = str(msg)
    if len(string) > length:
        return string[:length - 3] + " .."
    else:
        return string.ljust(length)
