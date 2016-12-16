#!/usr/bin/env python

import os.path, errno
import MDAnalysis as mda
import hop.graph, hop.walker
from hop.utilities import unlink_f, mkdir_p

import logging
logger = logging.getLogger('MDAnalysis.app')

def generate_walkertraj_locally(
