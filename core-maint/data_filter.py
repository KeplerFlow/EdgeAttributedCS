import random

def read_and_validate_edges(file_path, n):
    edges = []
    debug_set = set()
    
    with open(file_path, 'r') as file:
        first_line = file.readline()
        n_expected, m_expected = map(int, first_line.split())
        if n_expected != n:
            raise ValueError("Vertex count does not match the expected values.")
        
        for line in file:
            v1, v2 = map(int, line.split())
            # 忽略不符合要求的边，而不是抛出错误
            if v1 < 0 or v2 < 0 or v1 >= n or v2 >= n or v1 == v2:
                continue
            if v1 > v2:
                v1, v2 = v2, v1
            edge = (v1, v2)
            if edge in debug_set:
                continue
            debug_set.add(edge)
            edges.append(edge)
    
    random.shuffle(edges)
    return edges

def write_edges(edges, output_path):
    with open(output_path, 'w') as file:
        for v1, v2 in edges:
            file.write(f"{v1} {v2}\n")

if __name__ == "__main__":
    input_path = "/home/asc23/lcy/graph/Final/Dataset/email-Enron.txt"
    output_path = "/home/asc23/lcy/graph/Final/Dataset/email-Enron-coremain.txt"
    n = 36692  # You need to set this to the actual number of vertices
    
    try:
        edges = read_and_validate_edges(input_path, n)
        write_edges(edges, output_path)
        print("Edges processed and saved successfully.")
    except Exception as e:
        print(f"Error: {e}")
