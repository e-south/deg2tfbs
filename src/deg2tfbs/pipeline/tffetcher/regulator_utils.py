"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tffetcher/regulator_utils.py

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

def is_master_regulator(regulator: str) -> bool:
    master_regulators = {
        'flhd', 'flhc', 'fis', 'fhla', 'dksa', 'cysb', 'csgd', 'crp', 'cra', 'cpxr',
        'barr', 'argr', 'arca', 'srsr', 'spf', 'soxs', 'slya', 'ryhb', 'rob', 'rcsb',
        'purr', 'pohp', 'phob', 'pdhr', 'oxyr', 'ompr', 'nsrr', 'nhar', 'narp', 'narl',
        'nagc', 'nac', 'mode', 'mara', 'lrp', 'lexa', 'iscr', 'ihfb', 'ihfa', 'hns',
        'glng', 'glar', 'gcvb', 'gadx', 'gade', 'fur', 'fnr',
    }
    return regulator.lower().strip() in master_regulators

def is_sigma_factor(regulator: str) -> bool:
    sigma_factors = {
        'flia', 'feci', 'rpod', 'rpos', 'rpon', 'rpoh', 'rpoe'
    }
    return regulator.lower().strip() in sigma_factors

def is_nucleoid_regulator(regulator: str, nuc_reg_set: set) -> bool:
    reg_lower = regulator.lower().strip()
    return reg_lower in nuc_reg_set or any(reg_lower.startswith(nuc.lower()) for nuc in nuc_reg_set)
