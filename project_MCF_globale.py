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

    tipi = np.array([0], dtype = np.int8)
    energie = np.array([E0_MeV], dtype = np.float32)

    p_brem = 1 - math.exp(-step_s)
    p_pair = 1 - math.exp(-7 * step_s / 9)

    t_max = math.log(E0_MeV / E_el) / math.log(2)

    t = 0.0
    profile_t = []
    profile_n = []

    while t < t_det and len(tipi)>0:
        t += step_s

        if salva_sciame:
            profile_t.append(t)
            profile_n.append(len(tipi))

        mask_e = (tipi == 1)
        energie[mask_e] = dE_X0_MeV*step_s

        attivo =  energie >0
        tipi = tipi[attivo]
        energie = energie[attivo]

        if len(tipi) == 0:
          break

        mask_e = (tipi == 1) & (energie > E_el)
        mask_g = (tipi == 0) & (energie > soglia_prod)

        r = rng.random(len(tipi))

        do_brem = mask_e & (r<p_brem)
        do_pair = mask_g & (r<p_pair)

        new_tipi = []
        new_energie = []

        if np.any(do_brem):
          E2 = energie[do_brem]*0.5
          new_tipi.appen(no.ones(len(E2, dtype = np.int8)))
          new_tipi.appen(no.ones(len(E2, dtype = np.int8)))
          new_energie.appen(E2)
          new_energie.appen(E2)

        if np.any(do_pair):
          E2 = energie[do_pair]*0.5
          new_tipi.appen(no.ones(len(E2, dtype = np.int8)))
          new_tipi.appen(no.ones(len(E2, dtype = np.int8)))
          new_energie.appen(E2)
          new_energie.appen(E2)

        keep = (do_brem | do_pair)
        if np.any(keep):
            new_tipi.append(tipi[keep])
            new_energie.append(energie[keep])
        
        tipi = np.concatenate(new_tipi)
        energie = np.concatenate(new_energie)

        if not np.any((tipi==1)&(energie>E_el)) and 
             \not np.any ((tipi==0)&(energie>soglia_prod)):
               break
               
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
    plt.figure(figsize=(14,8))

    colors = {0: "red", 20: "green", 40: "blue"}

    for theta, hits in results.items():
        hits_pos = hits[hits > 0

        if len(hits_pos) < 10:
            continue
        xmin = hits.pos.min()
        xmax = hits.pos.max()

        if xmin <= 0:
          xmin = np.min(hits_pos[hits_pos>0])

        nbins = 30
        edges = np.logspace(np.log10(xmin), np.log10(xmax), nbins + 1)

        hits, edges = np.histogram(hits_pos, bins = edges, density = True)
        centres = np.sqrt(edges[:-1]*edges[1:])

        popt, _ = curve_fit(
            lognormal,
            centers,
            hist,
            p0=[np.log(np.mean(hits_pos)), 0.5],
            maxfev=20000
            )
        mu, sigma = popt

        x = np.logspace(np.log10(xmin),np.log10(xmax), 300)
        y = lognormal(x, mu_fit, sigma_fit)

        plt.hist(hits_pos, bins=edges, density=True, alpha=0.25,
                 color=colors[theta], label=f"istogramma {theta}°")

        plt.plot(x, lognormal(x,mu.sigma), color=colors[theta],
                 label=f"{theta}°fit")
    plt.xscale('log')
    plt.yscale('log')
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
    energies = np.logspace(0, 2, 5)
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

    E0_TeV = float(input("Energia iniziale (TeV), tra 1 e 100: "))
    step_s = float(input("Passo 0<s≤1: "))
    N_eventi = int(input("Numero eventi: "))

    E0_MeV = E0_TeV * 1e6

    print("\nProfilo sciame per 0°, 20°, 40")
    angoli = [0, 20, 40]
    colori = {0:"red", 20: "green", 40: "blue"}
    plt.figure(figsize=(8,5))

    step_profile = 0.01

    for theta in angoli:
      n, t_list, n_list = simulate_shower(E0_MeV, step_s_profile, theta, dE_X0, rng, salva_sciame = True)
      plt.plot(t_list, n_list, color = colori[theta], label = f"{theta}°")
      
    plt.xlabel("t (X0)")
    plt.ylabel("particelle")
    plt.title("profili dello sciame a diversi angoli")
    plt.grid(True)
    plt.legend()
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
