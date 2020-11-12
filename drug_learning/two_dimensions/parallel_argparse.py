def parse_args(parser):

    parser.add_argument("-r", "-ray",
                        dest = "ray",
                        action ='store_true',
                        help = "Parallelize with ray.")

    parser.add_argument('-ncpus', "--ncpus",
                        dest = "n_cpus",
                        action ='store',
                        default = None,
                        type = int,
                        help = "Number of cpus.")
