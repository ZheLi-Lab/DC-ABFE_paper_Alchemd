from alchemd.utils.file_parser import InputParser
from alchemd.utils.run import RunAlchemdSimulation
from optparse import OptionParser
import sys

class OptParser:  
    def __init__(self, fake_args=None):
        parser = OptionParser()
        parser.add_option('-p', '--topology', dest='topology', 
                         help='The filename of the toplogy file.', 
                         default='protein.prmtop')
        parser.add_option('-c', '--coordinate', dest='coordinate', 
                         help='The filename of the coordinate file.', 
                         default='protein.rst7')
        parser.add_option('-i', '--input', dest='input', 
                         help='The filename of the input file.', 
                         default='input.txt')
        if fake_args:
            self.option, self.args = parser.parse_args(fake_args)
        else:
            self.option, self.args = parser.parse_args()

def main():
    opts = OptParser()
    complex_coor = opts.option.coordinate
    complex_topo = opts.option.topology
    input_file = InputParser(opts.option.input)
    run = RunAlchemdSimulation(input_file, complex_coor, complex_topo)
    run.run()
    return 0

if __name__ == '__main__':
    sys.exit(main())
