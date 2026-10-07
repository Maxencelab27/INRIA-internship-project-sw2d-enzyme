import matplotlib.pyplot as plt
import numpy as np

L = 1000
N = 100
dx = L / N

cell_ids = np.arange(1, N + 1)
x_centers = (cell_ids - 0.5) * dx

h0 = np.where(x_centers < 500, 10.0, 1.0)
phi = np.where(x_centers < 500, 1.0, 0.1)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8.5), sharex=False)

def draw_walls(ax, is_meters=False):
    if is_meters:
        left_wall = 0.0
        right_wall = L
        text_left_x = 10.0
        text_right_x = L - 10.0
    else:
        left_wall = 0.5
        right_wall = N + 0.5
        text_left_x = 1.3
        text_right_x = N - 0.8

    ax.axvline(x=left_wall, color='darkred', linestyle='-', linewidth=3)
    ax.text(text_left_x, ax.get_ylim()[1]*0.5, 'Paroi\n(q = 0)', 
            color='darkred', fontsize=9, fontweight='bold', ha='left', va='center',
            bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="darkred", alpha=0.8))
    ax.axvline(x=right_wall, color='darkred', linestyle='-', linewidth=3)
    ax.text(text_right_x, ax.get_ylim()[1]*0.5, 'Paroi\n(q = 0)', 
            color='darkred', fontsize=9, fontweight='bold', ha='right', va='center',
            bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="darkred", alpha=0.8))

ax1.step(x_centers, h0, where='mid', color='dodgerblue', linewidth=2.5, label="Hauteur d'eau h0")
ax1.fill_between(x_centers, h0, step="mid", color='dodgerblue', alpha=0.15)

ax1.axvline(x=500.0, color='orange', linestyle='--', linewidth=2)
ax1.text(490.0, 5.5, 'Rupture de barrage\n(x = 500 m)', color='orange', ha='right', va='center', fontweight='bold')

ax1.set_ylabel("Hauteur d'eau $h_0$ (m)", fontsize=11, fontweight='bold')
ax1.set_ylim(0, 12)
ax1.set_xlim(0.0, L)
ax1.grid(True, linestyle=':', alpha=0.6)

ax1.set_title("État initial ($u_0 = 0$ m/s)", fontsize=13, fontweight='bold', color='black', loc='center', pad=12)

ax1.set_xlabel('Position physique x (m)', fontsize=11, fontweight='bold', color='dimgray', loc='left')

draw_walls(ax1, is_meters=True)

ax2.step(cell_ids, phi, where='mid', color='forestgreen', linewidth=2.5, label="Porosité cible phi*")
ax2.fill_between(cell_ids, phi, step="mid", color='forestgreen', alpha=0.15)

ax2.axvline(x=50.5, color='orange', linestyle='--', linewidth=2)
ax2.text(52.0, 0.55, 'Interface de porosité\n(Moitié poreuse à droite)', color='orange', ha='left', va='center', fontweight='bold')

ax2.set_ylabel('Porosité cible $\phi^*$', fontsize=11, fontweight='bold')
ax2.set_ylim(0, 1.2)
ax2.set_xlim(0.5, N + 0.5)
ax2.grid(True, linestyle=':', alpha=0.6)

ax2.set_title("Carte de porosité", fontsize=13, fontweight='bold', color='black', loc='center', pad=12)

ax2.set_xlabel('ID de la cellule (Discrétisation N = 100)', fontsize=11, fontweight='bold', color='dimgray', loc='left')

draw_walls(ax2, is_meters=False)

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)
plt.show()
