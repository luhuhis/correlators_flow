#!/usr/local/bin/python3.7m
import lib_process_data as lpd
import numpy


def main():
    parser, requiredNamed = lpd.get_parser()

    parser.add_argument('--PathOutputFolder', help='the path of output folder like /a/b', type=str)

    # identify the run
    requiredNamed.add_argument('--nsamples', help='number of bootstrap samples', type=int, required=True)

    args = parser.parse_args()

    # read in the normalized correlator
    inputfolder = "../"+lpd.get_merged_data_path(args.qcdtype, args.corr, "")+"/spf/"
    outputfolder = inputfolder if not args.PathOutputFolder else args.PathOutputFolder

    labels =      ("2_s2_alpha_a_4",    "2_s2_alpha_a_5",    "2_s1_alpha_a_4",   "2_s1_alpha_a_5",     "2_s2_beta_a_4",    "2_s2_beta_a_5",    "2_s1_beta_a_4",    "2_s1_beta_a_5")
    labels_plot = (r'$s_2 \alpha a 4$', r'$s_2 \alpha a 5$', r'$s_1 \alpha a 4$', r'$s_1 \alpha a 5$', r'$s_2 \beta a 4$', r'$s_2 \beta a 5$', r'$s_1 \beta a 4$', r'$s_1 \beta a 5$')
    kappas = []
    for label in labels:
        try:
            kappas.append(numpy.genfromtxt(inputfolder+"fitparams_chisqdof_"+label+"_"+str(args.nsamples)+".dat", skip_footer=1)[0:2, 0])
        except:
            print(label, " not found")
            kappas.append([numpy.nan, numpy.nan])
    print(kappas)

    pos = 0.3+numpy.asarray((0.45, 0.85, 1.25, 1.65, 2.65, 3.05, 3.45, 3.85))
    colors = ('red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue')

    xlabel = r'$\kappa / T^3$'

    fig, ax, plots = lpd.create_figure(xlabel=xlabel, xlabelpos=(0.5, -0.1), UseTex=False) #ylims=(0, 5), xlims=(0,5),

    for i, label in enumerate(labels):
        ax.errorbar(kappas[i][0], pos[i], xerr=kappas[i][1], color=colors[i], fmt='x-', fillstyle='none', markersize=5, mew=0.25, lw=0.8, elinewidth=0.5, capsize=1.2)

    ax.set_yticks(pos)
    ax.set_yticklabels(labels_plot)

    fig.savefig(outputfolder+"/kappas_"+args.corr+str(args.nsamples)+".pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
