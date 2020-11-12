def parse_args(parser):

    parser.add_argument("-s", "-split",
                        dest = "split",
                        action ='store_true',
                        help = "Split the input sdf files into chunks.")

    parser.add_argument('-ch', "--chunk",
                        dest = "n_chunks",
                        action ='store',
                        default = 1000,
                        type = int,
                        help = "Number of molecules in each chunk.")
