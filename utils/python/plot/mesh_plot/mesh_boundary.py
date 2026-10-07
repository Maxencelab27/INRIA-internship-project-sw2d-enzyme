import matplotlib.pyplot as plt
import sys
from collections import defaultdict

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


def find_boundary_edges(elements):
    """Retourne un ensemble d'arêtes appartenant à un seul élément."""
    edge_count = defaultdict(int)

    for elem in elements:
        # Liste des arêtes d’un triangle ou quadrilatère
        for i in range(len(elem)):
            a = elem[i]
            b = elem[(i+1) % len(elem)]
            edge = tuple(sorted((a, b)))  # tri pour éviter (n1,n2)/(n2,n1)
            edge_count[edge] += 1

    # Une arête est frontière si elle apparaît une seule fois
    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]

    return boundary_edges


def plot_boundary_edges(nodes, boundary_edges):
    fig, ax = plt.subplots(figsize=(10, 8))

    for n1, n2 in boundary_edges:
        (x1, y1) = nodes[n1]
        (x2, y2) = nodes[n2]
        ax.plot([x1, x2], [y1, y2], "-r", linewidth=1.0)  # arêtes en rouge

    ax.set_aspect('equal', 'box')
    ax.set_title("Arêtes de frontière du maillage")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    plt.show()


# -----------------------------
# Exemple d'utilisation
# -----------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage : python mesh_boundary.py file.2dm")
        sys.exit(1)

    file_path = sys.argv[1]
    nodes, elements = read_2dm(file_path)

    boundary_edges = find_boundary_edges(elements)
    print(f"{len(boundary_edges)} arêtes de frontière détectées")

    plot_boundary_edges(nodes, boundary_edges)
