# -*- coding: utf-8 -*-
"""
Created on 1 May 2018
@author: Dylan Jones

This module contains helper-methods for displaying results and figures as latex code

"""


def subplot_figures(names, subcaptions=None, caption=None, label=None, loc=None, intend=None):
    subcaptions = [""]*len(names) if subcaptions is None else subcaptions
    caption = "" if caption is None else caption
    label = "" if label is None else label
    loc = "h" if loc is None else loc
    intend = " "*4 if intend is None else " "*intend

    latex_str = r"\begin{figure}[" + loc + "]\n"
    latex_str += intend + r"\centering" + "\n"
    i = 0
    for cap, name in zip(subcaptions, names):
        latex_str += intend + r"\subfloat[" + cap + r"]{{\includegraphics[width=0.45\textwidth]{" + name + "}}} \n"
        if i == len(names) - 1:
            pass
        elif i % 2 == 0:
            latex_str += intend + r"\quad" + "\n"
        else:
            latex_str += intend + r"\\" + "\n"
        i += 1
    latex_str += intend + r"\caption{" + caption + "} \n"
    latex_str += intend + r"\label{" + label + "} \n"
    latex_str += r"\end{figure}"
    return latex_str


def graph(name, small=False, loc=None, caption=None):
    caption = "Caption" if caption is None else caption
    loc = "" if loc is None else loc
    label = "fig: " + name.replace("_", " ")

    latex_str = "\smallgraph" if small else "\graph"
    latex_str += "{" + loc + "}"
    latex_str += "{" + name + "}"
    latex_str += "{" + caption + "}"
    latex_str += "{\label{" + label + "}}"
    return latex_str


def table(rows, headers=None, locs=None, intend=None):
    locs = "c"*len(rows[0]) if locs is None else locs
    intend = " "*4 if intend is None else " "*intend

    latex_str = r"\begin{center}" + "\n"
    latex_str += intend + r"\begin{tabular}{" + locs + "} \n"
    if headers:
        latex_str += 2*intend + " & ".join(headers) + r"\\" + "\n" + 2*intend + r"\hline" + "\n"
    for i, row in enumerate(rows):
        str_row = " & ".join(row).replace(".", ",")
        latex_str += 2*intend + str_row
        if i < len(rows)-1:
            latex_str += r"\\"
        latex_str += "\n"

    latex_str += intend + r"\end{tabular}" + "\n" + r"\end{center}"
    return latex_str


def equation(name, func_str, label=None, intend=None):
    body = "equation*" if label is None else "equation"
    intend = " "*4 if intend is None else " "*intend
    lines = [name + "=" + func_str]
    if label:
        lines.append(r"\label{" + label + "}")
    return _insert_into_body(body, lines, intend)


def _insert_into_body(body, lines, intend):
    latex_str = r"\begin{" + body + "}\n"
    lines = [lines] if isinstance(lines, str) else lines
    for line in lines:
        latex_str += intend + line + "\n"
    latex_str += r"\end{" + body + "}"
    return latex_str
