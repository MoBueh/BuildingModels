import numpy as np
import matplotlib.pyplot as plt

def abklingfunktion(t, A0, tau):
    return A0 * np.exp(-t/tau)

# Zeitpunkte
t = np.linspace(0, 10, 100)

# Parameter der Abklingfunktion
A0 = 1.0   # Amplitude zu t=0
tau = 1.2  # Abklingzeit

# Berechnung der abklingenden Werte
A = abklingfunktion(t, A0, tau)

# Plot der Abklingfunktion


fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))

# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.plot(t, A)
# labels

# ax1.yaxis.set_major_formatter('{:.1f}'.format)
ax1.set_yticks(np.arange(-0.5,1.01,0.1), minor = True)
plt.plot([0,6], [-0.2, -0.2], c = 'tab:green')
# ax1.set_yticks(np.arange(-0.5,1.01,0.5), ['', 'Rücklauftemperatur', '', 'Vorlauftemperatur'])
ax1.set_yticks(np.arange(-0.5,1.01,0.5), ['']*4)
ax1.set_xlim(0,6)
ax1.set_xticks(range(0,7,1), ['']*7)
ax1.set_xlabel('Rohrlänge in m', fontsize = 10)
ax1.set_ylabel('Temperatur in °C', fontsize = 10)
ax1.text(0.15, 0.97, 'Vorlauftemperatur', fontsize = 8)
ax1.text(0.15, -0.03, 'Rücklauftemperatur', fontsize = 8)
ax1.text(0.15, -0.28, 'Raumtemperatur', fontsize = 8)
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.set_title(r'$\bf{Temperaturverlauf\ innerhalb\ des\ TABS}$')

plt.show()

