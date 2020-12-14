
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.
from typing import List, Union
from warnings import warn

import copy
from .utils import rev_dir

class ComponentEnumerator:
    def __init__(self, name:str,enumerator_type=""):
        """a component enumerator's job is to create new components in a process similar to mechanisms."""
        self.name = name
        self.enumerator_type = enumerator_type
    def update_components(self,component=None) -> List:
        """this will create new components based on the input component somehow
        The child class should implement this
        :return: empty list"""
        warn(f"Default update_components called for ComponentEnumerator = {self.name}.")
        return []
    def __repr__(self):
        return self.name