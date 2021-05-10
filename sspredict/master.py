
import argparse
from argparse import RawTextHelpFormatter
try:
    from sspredict.make_prediction.input_format import input_format 
    from sspredict.make_prediction.make_prediction import make_prediction
except:
    from make_prediction.input_format import input_format
    from make_prediction.make_prediction import make_prediction
def main():
    parser = argparse.ArgumentParser(description='Predict Solid Solution Strengthening.',
                                        formatter_class=RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-f', 
                        '--filein', 
                        required=False,
                        type=str,
                        help='Input JSON file containing necessary material data.\n'
                             'If not sure, use --format tag to check first.')
    group.add_argument('--format',
                        type=str,
                        required=False,
                        help='Use this tag to check general input JSON format for different models.\n'
                             'Available input formats: \n'
                             'FCC_BCC_Edge_Ternary\n'
                             'FCC_BCC_Edge_Composition_Temperature\n'
                             'BCC_Screw_Curtin_Ternary\n'
                             'BCC_Screw_Curtin_Composition_Temperature\n'
                             'BCC_Screw_Suzuki_Temperature\n'
                             'BCC_Screw_Suzuki_Ternary\n')
    args = parser.parse_args()


    if args.format:
        input_format(args.format).print_sample_input()
    
    if args.filein:
        model_predict = make_prediction(args.filein)
        model_predict.predict()
        model_predict.writeoutput()

if __name__ == '__main__':
    main()