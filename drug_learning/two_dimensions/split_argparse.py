def parse_args(parser):
    subparsers = parser.add_subparsers(dest = "split")

    parser_split = subparsers.add_parser("split", help = "Split the input sdf files into chunks" )

    parser_split.add_argument('-ch', "--chunk",
                        dest = "n_chunks",
                        action ='store',
                        default = 1000,
                        type = int,
                        help = "Number of molecules in each chunk.")

    parser_split.add_argument('-nsp', "--nworkers-sp",
                        dest = "nw_sp",
                        action ='store',
                        default = 1,
                        type = int,
                        help = "Number of workers to parallelize the split into chunks.")
