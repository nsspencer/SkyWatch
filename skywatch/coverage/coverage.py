from skywatch.coverage.coverage_definition import CoverageDefinition
from astropy.time import Time
import concurrent.futures

class Coverage:
    def __init__(self) -> None:
        self.coverage_definitions = list()
    
    def add_coverage_definition(self, coverage_definition: CoverageDefinition) -> 'Coverage':
        if not isinstance(coverage_definition, CoverageDefinition):
            raise TypeError("Coverage definition must be of type CoverageDefinition.")
        
        self.coverage_definitions.append(coverage_definition)
        return self

    def calculate(self, time: Time, )