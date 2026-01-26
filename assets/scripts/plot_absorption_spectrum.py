import subprocess
import numpy as np
import matplotlib.pyplot as plt
import plotext as pltxt
from scipy.signal.windows import gaussian
from scipy.ndimage import convolve1d

def parse_orca_spectrum(outname, grep_velocity=False):
    """
    Parse ORCA output file to extract absorption spectrum data using grep.
    
    Parameters:
    -----------
    outname : str
        Path to the ORCA output file
    
    Returns:
    --------
    wavelengths : list
        Wavelengths in nm
    oscillator_strengths : list
        Oscillator strengths (fosc)
    """

    look_for = 'VELOCITY' if grep_velocity else 'ELECTRIC'
    # cmd = f"sed -n '/ABSORPTION SPECTRUM VIA TRANSITION {look_for} DIPOLE MOMENTS/,/^$/p' {outname} | sed '$,/.*/p'"
    cmd = f"sed -n \'/ABSORPTION SPECTRUM VIA TRANSITION {look_for} DIPOLE MOMENTS/,/^$/p\' {outname} | awk \'BEGIN {{RS=\"\"; FS=\"\\n\"}} {{last=$0}} END {{print last}}\'"
    output = subprocess.getoutput(cmd)
    
    if not output:
        raise ValueError("Could not find absorption spectrum data in the file")
    
    names = []
    wavelengths = []
    oscillator_strengths = []
    
    # Parse each line of data
    for line in output.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 5 and '->' in line:
            try:
                names.append(" ".join(parts[0:3]))
                wavelength = float(parts[5])
                fosc = float(parts[6])
                wavelengths.append(wavelength)
                oscillator_strengths.append(fosc)
            except (ValueError, IndexError):
                continue
    
    if not wavelengths:
        raise ValueError("No valid transitions found in the file")
    
    return names, wavelengths, oscillator_strengths

def broaden_spectrum(wavelengths, oscillator_strengths, wl_range=(150, 300), 
                     num_points=1000, fwhm=10):
    """
    Create a broadened absorption spectrum using Gaussian functions.
    
    Parameters:
    -----------
    wavelengths : list
        Transition wavelengths in nm
    oscillator_strengths : list
        Oscillator strengths
    wl_range : tuple
        Wavelength range for the spectrum (min, max) in nm
    num_points : int
        Number of points in the spectrum
    fwhm : float
        Full width at half maximum for Gaussian broadening in nm
    
    Returns:
    --------
    wl_axis : array
        Wavelength axis
    spectrum : array
        Broadened absorption spectrum
    """
    wl_axis = np.linspace(wl_range[0], wl_range[1], num_points)
    spectrum = np.zeros(num_points)
    
    # Convert FWHM to standard deviation
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    
    # Add each transition as a Gaussian directly
    for wl, fosc in zip(wavelengths, oscillator_strengths):
        # Create Gaussian centered at this wavelength
        gaussian_peak = fosc * np.exp(-((wl_axis - wl)**2) / (2 * sigma**2))
        spectrum += gaussian_peak
    
    return wl_axis, spectrum

def plot_absorption_spectrum(outname, fwhm=10, 
                            show_sticks=True, print_to_terminal=True, save_fig=True):
    """
    Plot absorption spectrum from ORCA output file.
    
    Parameters:
    -----------
    outname : str
        Path to the ORCA output file
    wl_range : tuple
        Wavelength range for plotting (min, max) in nm
    fwhm : float
        Full width at half maximum for Gaussian broadening in nm
    show_sticks : bool
        Whether to show stick spectrum
    save_fig : str, optional
        Filename to save the figure
    """
    # Parse data
    names, wavelengths, oscillator_strengths = parse_orca_spectrum(outname)
    _, _, oscillator_strengths_velocities = parse_orca_spectrum(outname, grep_velocity=True)
    
    wl_range = (min(wavelengths)-20, max(wavelengths)+20)

    # Create broadened spectrum
    wl_axis, spectrum = broaden_spectrum(wavelengths, oscillator_strengths, 
                                        wl_range=wl_range, fwhm=fwhm)
    
    _, velocity_spectrum = broaden_spectrum(wavelengths, oscillator_strengths_velocities,
                                        wl_range=wl_range, fwhm=fwhm)

    norm_factor = max(spectrum)
    spectrum /= norm_factor
    velocity_spectrum /= max(velocity_spectrum)

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot broadened spectrum
    green = (65/255, 165/255, 165/255)
    lightgreen = (170/255, 225/255, 225/255)
    red = (232/255, 111/255, 136/255)
    ax.plot(wl_axis, spectrum, '-', color=green, linewidth=2, label='Line-broadened spectrum (electric dipole moments)')
    ax.plot(wl_axis, velocity_spectrum, '-', color=lightgreen, linewidth=2, label='Line-broadened spectrum (velocity dipole moments)')
    
    # Plot stick spectrum
    if show_sticks:
        for wl, fosc in zip(wavelengths, oscillator_strengths):
            if wl_range[0] <= wl <= wl_range[1]:
                ax.plot([wl, wl], [0, fosc/norm_factor], '-', color=red, linewidth=1.5, alpha=0.6)
        ax.plot([], [], '-', color=red, linewidth=1.5, alpha=0.6, label='Transitions (sticks, intensity from electric dipole moments)')
    
    # Formatting
    ax.set_xlabel('Wavelength (nm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Oscillator Strength', fontsize=12, fontweight='bold')
    ax.set_title('UV-Vis Absorption Spectrum', fontsize=14, fontweight='bold')
    ax.set_xlim(wl_range)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    
    if save_fig:
        svgname = outname.rstrip(".out") + "_spectrum.svg"
        plt.savefig(svgname, dpi=300, bbox_inches='tight')
        print(f"\n--> Saved figure to {svgname}")

    if print_to_terminal:
        pltxt.theme("pro")
        pltxt.plotsize(100,25)
        pltxt.cld()
        pltxt.plot(wl_axis, spectrum, color=37)
        pltxt.xlabel("Wavelength (nm)")
        pltxt.ylabel("Oscillator Strength")
        pltxt.show()

        print("\nWavelengths and Oscillator Strengths\n" \
              "from ELECTRIC dipole moments, normalized to 100%")
        print("------------------------------------------------")
        for name, wl, fosc in zip(names, wavelengths, oscillator_strengths):
            print(f"    {name:<15s}    {wl:.1f} nm     {fosc/max(oscillator_strengths)*100:>6.1f} %")

        print("\nWavelengths and Oscillator Strengths\n" \
              "from VELOCITY dipole moments, normalized to 100%")
        print("------------------------------------------------")
        for name, wl, fosc in zip(names, wavelengths, oscillator_strengths_velocities):
            print(f"    {name:<15s}    {wl:.1f} nm     {fosc/max(oscillator_strengths_velocities)*100:>6.1f} %")