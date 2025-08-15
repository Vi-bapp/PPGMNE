# 1. Inicialização
import gmsh
import numpy as np

# Inicializa o Gmsh.
gmsh.initialize()

# Opcional: define a visibilidade da interface gráfica (GUI).
gmsh.option.setNumber("General.Terminal", 1)

# ---

# 2. Carregamento e Extração da Malha
print("\n--- Carregando a malha do arquivo 'output.msh' ---")

# Carrega o arquivo de malha existente
try:
    gmsh.open("16x16_t.msh")
    print("Arquivo 'output.msh' carregado com sucesso.")
except Exception as e:
    print(f"Erro ao carregar o arquivo: {e}")
    gmsh.finalize()
    exit()

# Sincroniza o modelo (passo essencial antes de extrair a malha)
gmsh.model.geo.synchronize()

# ---

# 3. Extração e Geração dos Arquivos
print("\n--- Extração e Geração de Arquivos da Malha ---")

# --- Exporta as coordenadas (x, y) dos nós ---
nodeTags, coord, _ = gmsh.model.mesh.getNodes()
xyz_coordinates = coord.reshape(-1, 3)
xy_coordinates = xyz_coordinates[:, :2]
np.savetxt("coordenadas.dat", xy_coordinates, fmt='%f', delimiter=' ')
print("Arquivo 'coordenadas.dat' gerado com sucesso!")

# -----------------------------------------------

# --- Exporta a incidência nodal de cada elemento ---
# A dimensão do elemento (2D) e a tag da entidade geométrica (1)
elementTypes, elementTags, nodeTagsByElement = gmsh.model.mesh.getElements(2, 1)
if elementTypes:
    propriedades_do_elemento = gmsh.model.mesh.getElementProperties(elementTypes[0])
    num_nodes_per_elem = propriedades_do_elemento[3]
    connectivity_matrix = nodeTagsByElement[0].reshape(-1, num_nodes_per_elem)
    np.savetxt("incidencia_nodal.dat", connectivity_matrix, fmt='%i', delimiter=' ')
    print("Arquivo 'incidencia_nodal.dat' gerado com sucesso!")
else:
    print("Nenhum elemento 2D encontrado na superfície com tag 1.")

# -----------------------------------------------

# --- Exporta a numeração dos nós de contorno ---
boundary_entities = gmsh.model.getEntities(1)
boundary_node_tags = []
for entity in boundary_entities:
    dim = entity[0]
    tag = entity[1]
    nodes_on_entity = gmsh.model.mesh.getNodes(dim, tag)[0]
    boundary_node_tags.extend(nodes_on_entity)
boundary_node_tags = sorted(list(set(boundary_node_tags)))
boundary_node_array = np.array(boundary_node_tags)
np.savetxt("nos_de_contorno.dat", boundary_node_array, fmt='%i', delimiter=' ', newline=' ')
print("Arquivo 'nos_de_contorno.dat' gerado com sucesso!")

# ---

# 4. Finalização
gmsh.finalize()
