import os
import numpy as np
import glob
import ba_tools as ba
import matplotlib.pyplot as plt
import Klaas_Functions as func
import Klaas_Tools as Tools
import matplotlib.mlab as mlab

def default_pars():
    pars = {}
    pars['kT'] = 4.114  # pN nm
    pars['L0'] = 0.34  # nm / base pair
    pars['L_bp'] = 4753  # number of base pairs
    pars['P_nm'] = 50  # persistence length
    pars['S_pN'] = 1000  # stretch modulus
    pars['z0_nm'] = 0  # offset in nm / subunit
    pars['NRL'] = 168  # nucleosome repeat length
    pars['repeats'] = 16  # number of repeats
    pars['type'] = "Human"  # type of histone
    pars['NRL_str'] = str(pars['NRL'])+'x'+str(pars['repeats'])+'_'+pars['type']  # Nucleosome Repeat Length + #repeats
    pars['drift'] = []
    pars['save'] = True
    return pars

p = default_pars()

def read_fitfiles():
    fitfile_path = "C:\\Users\\brouw\\Desktop\\Data\\180604\\fitfiles\\"
    newpath = fitfile_path + r'\Figures'  # New path to save the figures
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    fitfiles = []
    os.chdir(fitfile_path)
    for file in glob.glob("*.fit"):
        fitfiles.append(file)

    print(fitfiles)
    steps, stacks = [], []  # used to save data (T-test)
    Steps, Stacks = [], []  # used to save data (Smoothening)
    F_Rup_up, Step_up, F_Rup_down, Step_down = [], [], [], []  # Rupture forces and corresponding jumps

    BT_Ruptures = np.empty((0, 3))  # Brower-Toland
    BT_Ruptures_Stacks = np.empty((0, 3))

    Fignum = 1  # Used for output line

    Filenames = []
    for Filenum, Filename in enumerate(fitfiles):
        print("Processing fitfile... " + str(Filename))

        f_pull, f_release, z_pull, z_release, z_fit_pull, transitions, time_pull, time_release, F, Z, T = ba.read_fitfiles(fitfile_path, Filename, p)
        f_wlc = np.logspace(np.log10(0.15), np.log10(int(np.max(f_pull))), 1000)


        LogFile = Tools.read_log(Filename[:-4] + '.log')  # loads the log file with the same name
        if LogFile: Pars = Tools.log_pars(LogFile)                                 #Reads in all the parameters from the logfile
        else:
            continue

        # Remove all datapoints that should not be fitted
        F_Selected = f_pull
        Z_Selected = z_pull * 1000
        T_Selected = time_pull

        if len(Z_Selected) < 10:
            print("<<<<<<<<<<<", Filename, '==> No data points left after filtering!>>>>>>>>>>>>')
            continue

        PossibleStates, ProbSum, Peak, States, AllStates, Statemask, NewStates, NewAllStates, NewStateMask = func.find_states_prob(
            F_Selected, Z_Selected, F, Z, Pars, MergeStates=True, Z_Cutoff=2)  # Finds States

        # Calculates stepsize
        Unwrapsteps = []
        Stacksteps = []
        for x in States:
            if x >= Pars['Fiber0_bp']:
                Unwrapsteps.append(x)
            else:
                Stacksteps.append(x)
        Stacksteps = func.state2step(Stacksteps)
        Unwrapsteps = func.state2step(Unwrapsteps)
        if len(Unwrapsteps) > 0: steps.extend(Unwrapsteps)
        if len(Stacksteps) > 0: stacks.extend(Stacksteps)

        # this plots the Force-Extension curve
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1, 2, 1)
        fig1.suptitle(Filename, y=.99)
        ax1.set_title(r'Extension-Force Curve')
        ax1.set_ylabel(r'Force (pN)')
        ax1.set_xlabel(r'Extension (nm)')
        ax1.scatter(Z, F, color='grey', lw=0.1, s=5)
        ax1.scatter(Z_Selected, F_Selected, color='black', lw=0.1, s=5)
        ax1.set_ylim([np.min(F_Selected) - 0.1 * np.max(F_Selected), np.max(F_Selected) + 0.1 * np.max(F_Selected)])
        ax1.set_xlim([np.min(Z_Selected) - 0.1 * np.max(Z_Selected), np.max(Z_Selected) + 0.1 * np.max(Z_Selected)])

        ax2 = fig1.add_subplot(1, 2, 2)
        ax2.set_title(r'Probability Landscape')
        ax2.set_xlabel(r'Contour Length (bp)')  # (nm) should be removed
        ax2.set_ylabel(r'Probability (AU)')
        ax2.plot(PossibleStates, ProbSum / np.sum(ProbSum), label='ProbSum')
        ax2.scatter(States, Peak / np.sum(ProbSum))

        # this plots the Timetrace
        fig2 = plt.figure()
        fig2.suptitle(Filename, y=.99)

        ax3 = fig2.add_subplot(1, 2, 1)
        ax3.set_title(r'Timetrace Curve')
        ax3.set_xlabel(r'Time (s)')
        ax3.set_ylabel(r'Extension (nm)')
        ax3.set_ylim([0, Pars['L_bp'] * Pars['DNAds_nm'] + 100])
        ax3.scatter(T, Z, color='grey', lw=0.1, s=5)

        ax4 = fig2.add_subplot(1, 2, 2, sharey=ax3)
        ax4.set_title(r'Probability Landscape')
        ax4.set_xlabel(r'Probability (AU)')
        ax4.plot(ProbSum, PossibleStates * Pars['DNAds_nm'])
        ax4.scatter(Peak, States * Pars['DNAds_nm'], color='blue')

        ax3.set_xlim([np.min(T_Selected) - 0.1 * np.max(T_Selected), np.max(T_Selected) + 0.1 * np.max(T_Selected)])
        ax3.set_ylim([np.min(Z_Selected) - 0.1 * np.max(Z_Selected), np.max(Z_Selected) + 0.1 * np.max(Z_Selected)])


        if len(States) < 1:
            print("<<<<<<<<<<<", Filename, '==> No States were found>>>>>>>>>>>>')
            continue

        ##############################################################################################
        ######## Begin Plotting Different States

        States = NewStates
        Statemask = NewStateMask
        AllStates = NewAllStates

        colors = [plt.cm.brg(each) for each in np.linspace(0, 1, len(States))]  # Color pattern for the states
        dX = 10  # Offset for text in plot

        # Calculate the rupture forces using a median filter
        a, b, c, d = func.RuptureForces(F_Selected, Z_Selected, T_Selected, States, Pars, ax1, ax3)
        F_Rup_up.extend(a)
        Step_up.extend(b)
        F_Rup_down.extend(c)
        Step_down.extend(d)

        # Brower-Toland analysis
        Rups = func.BrowerToland(F_Selected, Z_Selected, T_Selected, States, Pars, ax1, ax3)
        BT_Ruptures = np.append(BT_Ruptures, Rups, axis=0)

        # Brower-Toland analysis for stacking steps
        A = func.BrowerToland_Stacks(F_Selected, Z_Selected, T_Selected, States, Pars, ax1, ax3)
        BT_Ruptures_Stacks = np.append(BT_Ruptures_Stacks, A, axis=0)

        Sum = np.sum(Statemask, axis=1)
        ax1.scatter(Z_Selected[Sum == 0], F_Selected[Sum == 0], color='black',
                    s=20)  # Datapoint that do not belong to any state
        ax3.scatter(T_Selected[Sum == 0], Z_Selected[Sum == 0], color='black',
                    s=20)  # Datapoint that do not belong to any state

        # Plot the states and datapoints in the same color
        for j, col in zip(np.arange(len(colors)), colors):
            Mask = Statemask[:, j]
            Fit = AllStates[:, j]

            ax1.plot(Fit, F, alpha=0.9, linestyle=':', color=tuple(col))
            ax1.scatter(Z_Selected[Mask], F_Selected[Mask], color=tuple(col), s=20, alpha=.6)

            ax2.vlines(States[j], 0, np.max(Peak), linestyle=':', color=tuple(col))
            ax2.text(States[j], 0, int(States[j]), fontsize=10, horizontalalignment='center', verticalalignment='top',
                     rotation=90)

            ax3.plot(T, Fit, alpha=0.9, linestyle=':', color=tuple(col))
            ax3.scatter(T_Selected[Mask], Z_Selected[Mask], color=tuple(col), s=20, alpha=.6)

            ax4.hlines(States[j] * Pars['DNAds_nm'], 0, np.max(Peak), color=tuple(col), linestyle=':')
            ax4.text(0, States[j] * Pars['DNAds_nm'], int(States[j] * Pars['DNAds_nm']), fontsize=10,
                     verticalalignment='center', horizontalalignment='right')

        Unwrapsteps = []
        Stacksteps = []
        for x in NewStates:
            if x >= Pars['Fiber0_bp']:
                Unwrapsteps.append(x)
            else:
                Stacksteps.append(x)
        Stacksteps = np.diff(np.array(Stacksteps))
        Unwrapsteps = np.diff(np.array(Unwrapsteps))
        if len(Unwrapsteps) > 0: Steps.extend(Unwrapsteps)
        if len(Stacksteps) > 0: Stacks.extend(Stacksteps)
        ######################################################################################################################
        fig1.tight_layout()
        fig1.savefig(newpath + r'\\' + Filename[0:-4] + 'FoEx_all.png')

        fig2.tight_layout()
        fig2.savefig(newpath + r'\\' + Filename[0:-4] + 'Time_all.png')

        plt.close("all")

        Fignum += 2

        try:
            BT(BT_Ruptures, Pars, newpath, True)
        except:
            pass
        BT(BT_Ruptures_Stacks, Pars, newpath, False)

        # Plotting a histogram of the stepsizes
        fig3 = plt.figure()
        ax5 = fig3.add_subplot(1, 2, 1)
        ax6 = fig3.add_subplot(1, 2, 2)
        Range = [0, 400]
        Bins = 50
        n = ax5.hist(Steps, bins=Bins, range=Range, lw=0.5, zorder=1, color='blue', label='25 nm steps')[0]
        ax6.hist(Stacks, bins=int(Bins / 2), range=Range, lw=0.5, zorder=1, color='orange', label='Stacking transitions')

        # Fitting double gaussian over 25nm Steps
        try:
            Norm = Range[-1] / Bins
            D_Gaus = func.fit_2step_gauss(Steps)
            mu = D_Gaus[0]
            sigma = D_Gaus[1]
            x = np.linspace(mu - 3 * sigma, mu + 3 * sigma, 100)
            ax5.plot(x, mlab.normpdf(x, mu, sigma) * D_Gaus[2] * 2 * Norm, color='red', lw=4, zorder=10,
                     label='Gaussian fit')
            mu = 2 * mu
            x = np.linspace(mu - 3 * sigma, mu + 3 * sigma, 100)
            ax5.plot(x, mlab.normpdf(x, mu, sigma) * D_Gaus[3] * 2 * Norm, color='red', lw=4, zorder=10)
            ax5.text(Range[-1] - 100, np.max(n) - 0.1 * np.max(n), 'mean1:' + str(int(D_Gaus[0])),
                     verticalalignment='bottom')
            ax5.text(Range[-1] - 100, np.max(n) - 0.1 * np.max(n), 'mean2:' + str(int(2 * D_Gaus[0])),
                     verticalalignment='top')
        except:
            print('>>No 25 nm steps to fit gauss')

        ax5.set_xlabel('stepsize (bp)')
        ax5.set_ylabel('Count')
        ax5.set_title("Histogram stepsizes 25nm steps")
        ax5.legend(loc='best',
                   title='#Samples=' + str(len(Filenames)) + ', Binsize=' + str(int(np.max(Range) / Bins)) + 'bp/bin')
        ax6.set_xlabel('stepsize (bp)')
        ax6.set_ylabel('Count')
        ax6.set_title("Histogram stepsizes stacking steps")
        ax6.legend(loc='best', title='#Samples=' + str(len(Filenames)) + ', Binsize=' + str(
            int(np.max(Range) / int(Bins / 2))) + 'bp/bin')
        fig3.tight_layout()
        fig3.savefig(newpath + r'\\' + 'Hist.png')

        # plotting the rupture forces scatterplot
        fig4, ax7 = plt.subplots()
        ax7.scatter(F_Rup_up, Step_up, color='red', label='Jump to higher state')  # What should be the errors?
        ax7.scatter(F_Rup_down, Step_down, color='Green', label='Jump to lower state')  # What should be the errors?
        ax7.set_ylim(0, 400)
        ax7.set_xlabel('Rupture Forces (pN)')
        ax7.set_ylabel('Stepsize (bp)')
        ax7.set_title("Rupture forces versus stepsize")
        ax7.legend(loc='best')
        fig4.savefig(newpath + r'\\' + 'RF.png')

    return


def BT(BT_Ruptures, Pars, newpath, Steps=True):
    # Brower-Toland Analysis
    RFs = BT_Ruptures[:, 0]
    ln_dFdt_N = np.log(np.divide(BT_Ruptures[:, 2], BT_Ruptures[:, 1]))
    # Remove Ruptures at extensions larger than contour length (ln gets nan value)
    RFs = RFs[abs(ln_dFdt_N) < 10e6]
    ln_dFdt_N = ln_dFdt_N[abs(ln_dFdt_N) < 10e6]
    x = np.linspace(np.nanmin(ln_dFdt_N), np.nanmax(ln_dFdt_N), 10)
    a, a_err, b, b_err, d, D_err, K_d0, K_d0_err, Delta_G, Delta_G_err = func.dG_browertoland(ln_dFdt_N, RFs, Pars)

    # BowerToland plot
    fig, ax = plt.subplots()
    ax.plot(x, a * x + b, color='red', lw=2, label='Linear Fit')
    ax.plot(x, 1.3 * x + 19, color='green', lw=2, label='Result B-T')
    #    ax.plot(np.log(np.divide(A[:,2],A[:,1])), A[:,0], label='Data', color='red')
    ax.scatter(ln_dFdt_N, RFs, label='Data')
    ax.set_title("Brower-Toland analysis")
    Subtitle = "d = " + str(np.round(d, 1)) + "±" + str(np.round(D_err, 1)) + " nm"
    Subtitle = Subtitle + ", k_D(0) = {:.1e}".format(K_d0) + "±{:.1e}".format(K_d0_err) + " / sec"
    Subtitle = Subtitle + ", Delta G=" + str(Delta_G) + "±" + str(Delta_G_err) + " k_BT"
    fig.suptitle(Subtitle)
    ax.set_xlabel("ln[(dF/dt)/N (pN/s)]")
    ax.set_ylabel("Force (pN)")
    #    ax.set_ylim(5,40)
    #    ax.set_xlim(-4,2)
    ax.legend(loc='best',
              title='Slope:' + str(np.round(a, 1)) + '±' + str(np.round(a_err, 1)) + ', intersect:' + str(
                  np.round(b, 1)) + '±' + str(np.round(b_err, 1)))
    if Steps:
        fig.savefig(newpath + r'\\' + 'BT_Steps.png')
    else:
        fig.savefig(newpath + r'\\' + 'BT_Stacks.png')

    return


if __name__ == "__main__":
    steps_fitfiles()
