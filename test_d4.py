#!/usr/bin/env python3
"""Test script for D4 equivalence classes"""

from itertools import permutations

def apply_d4_operations(config):
    a, b, c, d = config
    return [
        (a, b, c, d),      # e: identity
        (d, a, b, c),      # r90: rotate 90°
        (c, d, a, b),      # r180: rotate 180°
        (b, c, d, a),      # r270: rotate 270°
        (d, c, b, a),      # sx: horizontal reflection
        (b, a, d, c),      # sy: vertical reflection
        (a, d, c, b),      # sd: main diagonal reflection
        (c, b, a, d),      # sd': anti-diagonal reflection
    ]

def canonical_form(config):
    return min(apply_d4_operations(config))

def find_d4_equivalence_classes(values):
    all_perms = list(permutations(values))
    canonical_forms = set()
    for perm in all_perms:
        canonical_forms.add(canonical_form(perm))
    return sorted(list(canonical_forms))

# Test with [1, 2, 3, 6]
test_values = [1, 2, 3, 6]
result = find_d4_equivalence_classes(test_values)
print(f"Input: {test_values}")
print(f"Number of equivalence classes: {len(result)}")
print(f"Equivalence classes:")
for i, config in enumerate(result, 1):
    print(f"  {i}. {list(config)}")
