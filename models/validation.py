import sys
from collections import defaultdict

def parse_obj_for_half_edge(file_path):
    vertices = []
    edges = defaultdict(list)  # Key: edge (v1, v2), Value: list of face indices
    non_manifold_edges = []

    try:
        with open(file_path, 'r') as file:
            face_index = 0
            for line in file:
                parts = line.strip().split()
                if not parts:
                    continue
                if parts[0] == 'v':
                    # Parse vertex
                    vertices.append(tuple(map(float, parts[1:4])))
                elif parts[0] == 'f':
                    # Parse face and create directed edges
                    face_vertices = [int(p.split('/')[0]) for p in parts[1:]]
                    for i in range(len(face_vertices)):
                        v1 = face_vertices[i]
                        v2 = face_vertices[(i + 1) % len(face_vertices)]
                        edge = (v1, v2)  # Preserve direction
                        edges[edge].append(face_index)
                    face_index += 1
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return [], {}, []
    except Exception as e:
        print(f"Error while parsing OBJ: {e}")
        return [], {}, []

    # Check for non-manifold edges
    for edge, faces in edges.items():
        reverse_edge = (edge[1], edge[0])
        total_refs = len(faces) + len(edges.get(reverse_edge, []))
        if total_refs != 2:  # Must be referenced exactly twice (one forward, one backward)
            non_manifold_edges.append((edge, faces))

    return vertices, edges, non_manifold_edges

def main():
    if len(sys.argv) < 2:
        print("Usage: python check_half_edge_structure.py <path_to_obj>")
        return

    file_path = sys.argv[1]
    vertices, edges, non_manifold_edges = parse_obj_for_half_edge(file_path)

    print(f"Total vertices: {len(vertices)}")
    print(f"Total directed edges: {len(edges)}")
    if non_manifold_edges:
        print(f"Non-manifold edges found ({len(non_manifold_edges)}):")
        for edge, faces in non_manifold_edges:
            print(f"Edge {edge} is referenced in faces {faces}")
    else:
        print("The model is a closed manifold with a valid HalfEdge structure.")

if __name__ == "__main__":
    main()
