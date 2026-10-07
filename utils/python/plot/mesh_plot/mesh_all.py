import matplotlib.pyplot as plt
import sys

def read_2dm(file_path):
    nodes = {}
    elements = []

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            # --- NODES ---
            if parts[0] == "ND":
                # Format : ND id x y z
                node_id = int(parts[1])
                x = float(parts[2])
                y = float(parts[3])
                nodes[node_id] = (x, y)

            # --- TRIANGLES (E3T) ---
            elif parts[0] == "E3T":
                # Format : E3T id n1 n2 n3 matid
                _, elem_id, n1, n2, n3, *_ = parts
                elements.append([int(n1), int(n2), int(n3)])

            # --- QUADS (E4Q) ---
            elif parts[0] == "E4Q":
                # Format : E4Q id n1 n2 n3 n4 matid
                _, elem_id, n1, n2, n3, n4, *_ = parts
                elements.append([int(n1), int(n2), int(n3), int(n4)])

    return nodes, elements


def plot_mesh(nodes, elements):
    fig, ax = plt.subplots(figsize=(10, 8))

    for elem in elements:
        # Récupère les points (x, y)
        pts = [nodes[n] for n in elem]
        # Boucle pour refermer le polygone
        pts.append(pts[0])

        xs, ys = zip(*pts)
        ax.plot(xs, ys, "-k", linewidth=0.5)

    ax.set_aspect('equal', 'box')
    ax.set_title("Maillage 2DM")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    plt.show()


# -----------------------------
# Exemple d'utilisation
# -----------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage : python mesh_all.py file.2dm")
        sys.exit(1)

    file_path = sys.argv[1]
    nodes, elements = read_2dm(file_path)
    plot_mesh(nodes, elements)
