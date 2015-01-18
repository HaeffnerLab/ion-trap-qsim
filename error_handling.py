#!/usr/bin/env python

class Dimensionerror(RuntimeError):

   def __init__(self, arg):

      self.args = arg + "Dimension of Hilbert Space must be greater than quantum state number!"

