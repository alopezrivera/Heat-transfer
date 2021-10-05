import matplotlib.pyplot as plt


def panel(*args):

    fs = 10  # sets the fontsize

    axes = [plt.subplot(int(f'22{i}')) for i in range(1, len(args)+1)]

    for i in range(len(axes)):
        axes[i].plot(args[i])

    plt.xlabel(r'$time \ [s]$', fontsize=fs)
    plt.ylabel(r'$pressure ratio \ [-]$', fontsize=fs)
    plt.grid(True)
    plt.show()
