import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# costanti fisiche
h_sciame = 20000
h_rilevatore = 4000
X0 = 700

strada = h_sciame - h_rilevatore
lunghezza_di_radiazione = strada / X0

E_el = 87.92
riposo_el = 0.511
soglia_prod = 2 * riposo_el
dE_X0 = 73.24


# simulazione sciame
def simulate_shower(E0_MeV, step_s=0.1, theta_deg=0.0, dE_X0_MeV=dE_X0,
                    rng=None, salva_sciame=False):

    if rng is None:
        rng = np.random.default_rng()

    theta_rad = np.deg2rad(theta_deg)
    t_det = lunghezza_di_radiazione / np.cos(theta_rad)

    particles_type = ['g']
    particles_E = [E0_MeV]

    p_brem = 1 - math.exp(-step_s)
    p_pair = 1 - math.exp(-7 * step_s / 9)

    t_max = math.log(E0_MeV / E_el) / math.log(2)

    t = 0.0
    profile_t = []
    profile_n = []

    while t < t_det and particles_type:
        t += step_s

        if salva_sciame:
            profile_t.append(t)
            profile_n.append(len(particles_type))

        new_type = []
        new_E = []

        attivo = (t <= t_max)

        for p_type, E in zip(particles_type, particles_E):

            if p_type == 'e':
                E -= dE_X0_MeV * step_s
                if E <= 0:
                    continue

                if attivo and E > E_el:
                    if rng.random() < p_brem:
                        E2 = E * 0.5
                        new_type += ['e', 'g']
                        new_E += [E2, E2]
                        continue

                new_type.append('e')
                new_E.append(E)

            else:  # fotone
                if E <= soglia_prod:
                    continue

                if attivo:
                    if rng.random() < p_pair:
                        E2 = E * 0.5
                        new_type += ['e', 'e']
                        new_E += [E2, E2]
                    else:
                        new_type.append('g')
                        new_E.append(E)
                else:
                    new_type.append('g')
                    new_E.append(E)

        particles_type = new_type
        particles_E = new_E

    if salva_sciame:
        return len(particles_type), profile_t, profile_n
    else:
        return len(particles_type)


# spettro E^-2
def sample_E_spectrum_Eminus2(N, Emin, Emax, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    u = rng.random(N)
    inv_Emin = 1 / Emin
    inv_Emax = 1 / Emax
    return 1 / (inv_Emin - u * (inv_Emin - inv_Emax))


# risposta rivelatore
def detector_response_flux(N_eventi=1000, step_s=0.1,
                           angles_deg=(0, 20, 40),
                           Emin_TeV=1, Emax_TeV=100,
                           dE_X0_MeV=dE_X0, rng=None):

    if rng is None:
        rng = np.random.default_rng()

    Emin_MeV = Emin_TeV * 1e6
    Emax_MeV = Emax_TeV * 1e6

    E_singola = sample_E_spectrum_Eminus2(N_eventi, Emin_MeV, Emax_MeV, rng)

    results = {}
    for theta in angles_deg:
        hits = []
        for E0 in E_singola:
            n = simulate_shower(E0, step_s, theta, dE_X0_MeV, rng)
            hits.append(n)
        results[theta] = np.array(hits)

    return results


# fit lognormale
def lognormal(x, mu, sigma):
    return 1.0 / (x * sigma * np.sqrt(2*np.pi)) * np.exp(
        - (np.log(x) - mu)**2 / (2 * sigma**2)
    )


def plot_lognormal_fits_with_hist(results):
    plt.figure(figsize=(12,7))

    colors = {0: "red", 20: "green", 40: "blue"}

    for theta, hits in results.items():
        hits_pos = hits[hits > 0]

        hist, bin_edges = np.histogram(hits_pos, bins=30, density=True)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        popt, pcov = curve_fit(
            lognormal,
            bin_centers,
            hist,
            p0=[np.log(np.mean(hits_pos)), 0.5]
        )
        mu_fit, sigma_fit = popt

        x = np.linspace(min(hits_pos), max(hits_pos), 300)
        y = lognormal(x, mu_fit, sigma_fit)

        plt.hist(hits_pos, bins=30, density=True, alpha=0.25,
                 color=colors[theta], label=f"istogramma {theta}°")

        plt.plot(x, y, color=colors[theta], linewidth=2,
                 label=f"fit {theta}° (mu={mu_fit:.2f}, sigma={sigma_fit:.2f})")

    plt.xlabel("hit")
    plt.ylabel("densità")
    plt.title("istogrammi + fit lognormali")
    plt.grid(True)
    plt.legend()
    plt.show()


# legge di potenza
def power_law(E, A, alpha):
    return A * E**alpha


def fit_power_law(step_s=0.1, theta=0, N=200):
    energies = np.logspace(0, 2, 6)
    means = []

    rng = np.random.default_rng()

    for E in energies:
        hits = []
        for _ in range(N):
            n = simulate_shower(E*1e6, step_s, theta, dE_X0, rng)
            hits.append(n)
        means.append(np.mean(hits))

    energies = np.array(energies)
    means = np.array(means)

    popt, pcov = curve_fit(power_law, energies, means)
    A_fit, alpha_fit = popt
    sigma_alpha = np.sqrt(pcov[1,1])

    E_plot = np.logspace(0, 2, 200)

    plt.figure(figsize=(8,5))
    plt.loglog(energies, means, 'o')
    plt.loglog(E_plot, power_law(E_plot, A_fit, alpha_fit),
               label=f"fit: alpha={alpha_fit:.2f}±{sigma_alpha:.2f}")
    plt.xlabel("energia [TeV]")
    plt.ylabel("hit medi")
    plt.title("legge di potenza")
    plt.legend()
    plt.grid(True)
    plt.show()


# fluttuazioni relative
def study_relative_fluctuations(step_s=0.1, theta=0, N=200):
    energies = np.logspace(0, 2, 6)
    rel_flucts = []

    rng = np.random.default_rng()

    for E in energies:
        hits = []
        for _ in range(N):
            n = simulate_shower(E*1e6, step_s, theta, dE_X0, rng)
            hits.append(n)
        hits = np.array(hits)
        rel_flucts.append(hits.std() / hits.mean())

    plt.figure(figsize=(8,5))
    plt.semilogx(energies, rel_flucts, 'o-')
    plt.xlabel("energia [TeV]")
    plt.ylabel("sigma/mu")
    plt.title("fluttuazioni relative")
    plt.grid(True)
    plt.show()


# main
def main():
    rng = np.random.default_rng()

    print("Simulazione sciami\n")

    E0_TeV = float(input("Energia iniziale (TeV): "))
    step_s = float(input("Passo 0≤s≤1: "))
    N_eventi = int(input("Numero eventi: "))

    E0_MeV = E0_TeV * 1e6

    print("\nProfilo sciame")
    n, t_list, n_list = simulate_shower(E0_MeV, step_s, 0, dE_X0, rng, salva_sciame=True)

    plt.figure(figsize=(8,5))
    plt.plot(t_list, n_list)
    plt.xlabel("t (X0)")
    plt.ylabel("particelle")
    plt.title("profilo sciame")
    plt.grid(True)
    plt.show()

    print("\nStudio statistico")
    results = detector_response_flux(
        N_eventi=N_eventi,
        step_s=step_s,
        angles_deg=(0, 20, 40),
        rng=rng
    )

    for theta, hits in results.items():
        print(f"\nAngolo {theta}°")
        print("media:", hits.mean())
        print("mediana:", np.median(hits))
        print("eventi con almeno 1 hit:", (hits > 0).mean())

    print("\nFit lognormale")
    plot_lognormal_fits_with_hist(results)

    print("\nFit legge potenza")
    fit_power_law(step_s=step_s)

    print("\nFluttuazioni relative")
    study_relative_fluctuations(step_s=step_s)


if __name__ == "__main__":
    main()
