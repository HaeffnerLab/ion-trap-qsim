#!/usr/bin/env python
from __future__ import division, absolute_import, print_function, unicode_literals

class Dimensionerror(RuntimeError):

   def __init__(self, *arg):

      self.args = ' '.join([arg, "Dimension of Hilbert Space must be greater than quantum state number!"])

