# 1. Inicialização
import gmsh
import numpy as np

# Inicializa o Gmsh.
# É o primeiro passo para usar a API.
gmsh.initialize()

# Opcional: define a visibilidade da interface gráfica (GUI).
gmsh.option.setNumber("General.Terminal", 1)
#gmsh.fltk.run() # Executa a interface gráfica para visualização da malha.

# ---

# 2. Criação da Malha
gmsh.model.add("my_model")
lc = 1.0 / 15.0 # Exemplo de tamanho característico da malha

# Adiciona os pontos da geometria
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(1, 0, 0, lc, 2)
gmsh.model.geo.addPoint(1, 1, 0, lc, 3)
gmsh.model.geo.addPoint(0, 1, 0, lc, 4)

# Conecta os pontos para formar as linhas
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

# Cria o loop de curvas e a superfície
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

# Sincroniza o modelo
gmsh.model.geo.synchronize()

# --- Aqui você define o tipo de elemento ---
# Ativa a recombinação para gerar malha de elementos quadriláteros.
# Comente ou remova esta linha para ter uma malha triangular.
gmsh.option.setNumber("Mesh.RecombineAll", 1)

# Gera a malha (neste ponto, ela é linear, seja triangular ou quadrilátera).
gmsh.model.mesh.generate(2)

# --- Aqui você define a ordem do elemento ---
# Converte a malha para segunda ordem (quadrática).
# Comente ou remova esta linha para ter uma malha linear.
#gmsh.model.mesh.setOrder(2)

# Opcional: Salva a malha para visualização
gmsh.write("output.msh")

# --- Abre a interface gráfica do Gmsh para visualização ---
# Remova o '#' para ativar a visualização da malha.
# Se você tiver problemas de visualização (erros X11),
# basta comentar esta linha novamente.
gmsh.fltk.run()

# ---

# 3. Extração e Geração dos Arquivos
print("\n--- Extração e Geração de Arquivos da Malha ---")

# --- Exporta as coordenadas (x, y) dos nós ---
# Extrai as coordenadas dos nós
nodeTags, coord, _ = gmsh.model.mesh.getNodes()

# Converte o array de coordenadas em uma matriz (N x 3) e seleciona as colunas X e Y
xyz_coordinates = coord.reshape(-1, 3)
xy_coordinates = xyz_coordinates[:, :2]

# Cria e salva o arquivo de coordenadas
np.savetxt("coordenadas.dat", xy_coordinates, fmt='%f', delimiter=' ')
print("Arquivo 'coordenadas.dat' gerado com sucesso!")

# -----------------------------------------------

# --- Exporta a incidência nodal de cada elemento ---
# Extrai a conectividade dos elementos
elementTypes, elementTags, nodeTagsByElement = gmsh.model.mesh.getElements(2, 1)

# Obtém o número de nós por elemento usando getElementProperties
propriedades_do_elemento = gmsh.model.mesh.getElementProperties(elementTypes[0])
num_nodes_per_elem = propriedades_do_elemento[3]

# Reorganiza o array de conectividade em uma matriz (N x num_nodes_per_elem).
connectivity_matrix = nodeTagsByElement[0].reshape(-1, num_nodes_per_elem)

# Cria e salva o arquivo de incidência nodal
np.savetxt("incidencia_nodal.dat", connectivity_matrix, fmt='%i', delimiter=' ')
print("Arquivo 'incidencia_nodal.dat' gerado com sucesso!")

# -----------------------------------------------

# --- Exporta a numeração dos nós de contorno ---
# Obtém as entidades de contorno (as linhas, que são entidades 1D) do modelo.
# O '1' significa dimensão 1D.
boundary_entities = gmsh.model.getEntities(1)

boundary_node_tags = []

# Itera sobre cada entidade de contorno (cada linha).
for entity in boundary_entities:
    dim = entity[0]
    tag = entity[1]
    
    # Obtém os nós da malha associados a esta entidade (linha).
    # Acessa os nós da entidade usando o método 'getNodes'
    nodes_on_entity = gmsh.model.mesh.getNodes(dim, tag)[0]
    boundary_node_tags.extend(nodes_on_entity)

# Remove duplicatas (nós compartilhados nas quinas) e ordena
boundary_node_tags = sorted(list(set(boundary_node_tags)))

# Cria um array NumPy a partir da lista
boundary_node_array = np.array(boundary_node_tags)

# Salva a numeração dos nós de contorno no arquivo
# 'delimiter' e 'newline' são definidos como um espaço em branco.
np.savetxt("nos_de_contorno.dat", boundary_node_array, fmt='%i', delimiter=' ', newline=' ')

print("Arquivo 'nos_de_contorno.dat' gerado com sucesso!")

# ---

# 4. Finalização
# Libera os recursos do Gmsh.
gmsh.finalize()