from matplotlib import pyplot as plt
import numpy as np


def plotPositions(txPosition, rxPosition, txNormal, rxNormal):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    # plot positions
    ax.plot(*txPosition, 'r*')
    ax.plot(*rxPosition, 'bx')

    # plot pointing vectors
    # looks like: ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)
    l = np.linalg.norm(txPosition - rxPosition)
    ax.quiver(*np.concatenate((txPosition, txNormal)),
              length=0.4 * l, normalize=False)
    # ax.quiver(*np.concatenate((rxPosition,rxNormal)), length=0.4*l, normalize=False)
    ax.quiver(*np.concatenate((txPosition, -rxNormal)),
              length=0.4 * l, normalize=False)

    # plot config
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 10])
    ax.set_zlim([-0.01, 10])
    plt.title('Satellite and RA Positions')
    return fig, ax


def plotAntennaPattern(antennaPatternFunction):
    N = 101
    az = np.linspace(-np.pi, np.pi, 101)
    el = np.linspace(-np.pi, np.pi, 101)
    AZ, EL = np.meshgrid(az, el)
    Z = antennaPatternFunction(AZ, EL)

    # plotting
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(AZ * (360 / (2 * np.pi)),
                           EL * (360 / (2 * np.pi)), Z, linewidth=0)
    plt.xlabel('Azimuth (Degrees)')
    plt.ylabel('Elevation (Degrees)')
    return fig, ax
