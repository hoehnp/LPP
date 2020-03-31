# movie of self-assembling micelles

from __future__ import absolute_import
from builtins import range

d = dump("dump.micelle")

s = svg(d)
s.acol([1,2,3,4],["blue","red","cyan","yellow"])
s.arad(list(range(4)),0.5)
s.zoom(1.5)

s.file = "micelle"

s.all()
