#!/usr/bin/env python
from __future__ import division, absolute_import, print_function, unicode_literals

class Dimensionerror(RuntimeError):

   def __init__(self, *arg):

      self.args = ' '.join([arg, "Dimension of Hilbert Space must be greater than quantum state number!"])



class Statetypeerror( Exception ):


        def __init__(self, *args):
                """
                Give error when the initialized state type is not what is supposed to be. Keyword arguments are ion_number and state_type
                For example, if the ion state type was set as pure and was initialized as density_operator.

                """

                self.args = "State of ion number %i must be a %s" % (args[0], args[1] )

        def __str__(self):
                return ''.join(self.args)


