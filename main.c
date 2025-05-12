#include <math.h>
#include <stdlib.h>  
#include <stdio.h>  
#include <string.h> 

typedef struct Vertex Vertex;
typedef struct HalfEdge HalfEdge;
typedef struct Face Face;   

struct Vertex {
    int id;
    int x, y;
    HalfEdge *incident_edge;  // Uma half-edge que parte deste vértice
};

struct HalfEdge {
    int id;
    Vertex *origin;           // Vértice de origem
    HalfEdge *twin;           // Half-edge oposta
    HalfEdge *next;           // Próxima half-edge na face
    HalfEdge *prev;           // Half-edge anterior na face (opcional, mas útil)
    Face *face;               // Face à esquerda dessa half-edge
};

struct Face {
    int id;
    HalfEdge *outer_component; // Uma half-edge que contorna a face
};

typedef struct {
    int u, v;
    HalfEdge *edge;
} EdgeRecord;


typedef struct {
    int x, y;
} Coord;

typedef struct {
    int *vertex_indices;  // Lista de índices de vértices que formam a face
    int n_vertices;       // Quantidade de vértices na face
} RawFace;

typedef struct {
    int u, v;
    int count;  // Número de ocorrências
} EdgeCount;


EdgeCount *edge_counts = NULL;
int ec_size = 0;
int ec_cap = 0;

EdgeRecord *edge_map = NULL;
int edge_map_size = 0;
int edge_map_capacity = 0;



// Função para alocar leitura segura de vértices
Coord *read_vertices(int n) {
    Coord *vertices = malloc(sizeof(Coord) * n);
    if (!vertices) {
        perror("Erro ao alocar vértices");
        exit(EXIT_FAILURE);
    }


    for (int i = 0; i < n; i++) {
        scanf("%d %d", &vertices[i].x, &vertices[i].y);
    }
    return vertices;
}

// Função para ler as faces com listas variáveis de vértices
RawFace *read_faces(int f) {
    RawFace *faces = malloc(f * sizeof(RawFace));
    if (!faces) {
        perror("Erro ao alocar faces");
        exit(EXIT_FAILURE);
    }

    char buffer[4096];
    for (int i = 0; i < f; ++i) {
        // Lê a linha inteira
        fgets(buffer, sizeof(buffer), stdin);  // pode haver \n pendente da leitura anterior
        while (buffer[0] == '\n') fgets(buffer, sizeof(buffer), stdin);  // pula linhas vazias

        // Conta quantos números tem
        int count = 0;
        for (char *p = buffer; *p; ++p) {
            if (*p == ' ') count++;
        }
        count++;  // número de vértices na face

        // Aloca e lê os índices
        faces[i].vertex_indices = malloc(count * sizeof(int));
        faces[i].n_vertices = count;

        char *p = buffer;
        for (int j = 0; j < count; ++j) {
            sscanf(p, "%d", &faces[i].vertex_indices[j]);
            while (*p && *p != ' ') p++;  // avança até espaço
            while (*p == ' ') p++;        // pula espaços
        }
    }
    return faces;
}

void add_or_increment_edge(int u, int v) {
    for (int i = 0; i < ec_size; ++i) {
        if (edge_counts[i].u == u && edge_counts[i].v == v) {
            edge_counts[i].count++;
            return;
        }
    }

    if (ec_size >= ec_cap) {
        ec_cap = ec_cap ? 2 * ec_cap : 64;
        edge_counts = realloc(edge_counts, ec_cap * sizeof(EdgeCount));
    }

    edge_counts[ec_size++] = (EdgeCount){ u, v, 1 };
}

const char *check_edge_counts(RawFace *faces, int n_faces) {
    for (int i = 0; i < n_faces; ++i) {
        int nv = faces[i].n_vertices;
        int *idx = faces[i].vertex_indices;
        for (int j = 0; j < nv; ++j) {
            int u = idx[j] - 1;
            int v = idx[(j + 1) % nv] - 1;
            add_or_increment_edge(u, v);
        }
    }

    for (int i = 0; i < ec_size; ++i) {
        int u = edge_counts[i].u;
        int v = edge_counts[i].v;

        // Conta quantas vezes a aresta (v → u) também aparece
        int total = edge_counts[i].count;
        for (int j = 0; j < ec_size; ++j) {
            if (edge_counts[j].u == v && edge_counts[j].v == u)
                total += edge_counts[j].count;
        }

        if (total < 2) return "aberta";
        if (total > 2) return "nao subdivisao planar";
    }

    return NULL; // OK
}

// Verifica se dois segmentos [p1, p2] e [q1, q2] se intersectam
int segments_intersect(Coord p1, Coord p2, Coord q1, Coord q2) {
    int det1 = (q2.x - q1.x) * (p1.y - q1.y) - (q2.y - q1.y) * (p1.x - q1.x);
    int det2 = (q2.x - q1.x) * (p2.y - q1.y) - (q2.y - q1.y) * (p2.x - q1.x);
    int det3 = (p2.x - p1.x) * (q1.y - p1.y) - (p2.y - p1.y) * (q1.x - p1.x);
    int det4 = (p2.x - p1.x) * (q2.y - p1.y) - (p2.y - p1.y) * (q2.x - p1.x);
    return (det1 * det2 < 0) && (det3 * det4 < 0);// estritamente dentro dos segmentos
}

const char *check_face_intersections(Coord *vertices, RawFace *faces, int n_faces) {
    for (int i = 0; i < n_faces; ++i) {
        // printf("Face: %d    ", i);
        RawFace *f1 = &faces[i];
        int *idx1 = f1->vertex_indices;
        int n1 = f1->n_vertices;

        for (int j = i + 1; j < n_faces; ++j) {
            // printf("Comparando com Face: %d\n", j);
            RawFace *f2 = &faces[j];
            int *idx2 = f2->vertex_indices;
            int n2 = f2->n_vertices;

            for (int a = 0; a < n1; ++a) {
                int a1 = idx1[a] - 1;
                int a2 = idx1[(a + 1) % n1] - 1;

                // printf("A1:%d A2:%d  ", a1, a2);

                for (int b = 0; b < n2; ++b) {
                    int b1 = idx2[b] - 1;
                    int b2 = idx2[(b + 1) % n2] - 1;
                    
                    // printf("Comparando com B1:%d B2:%d\n", b1, b2);


                    // Ignora segmentos que compartilham vértices
                    if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;

                    int intersect = segments_intersect(vertices[a1], vertices[a2], vertices[b1], vertices[b2]);
                    
                    if (intersect) {
                        // printf("Segmento 1:(%d,%d) (%d,%d)  Segmento 2:(%d,%d) (%d,%d)",  vertices[a1].x, vertices[a1].y, vertices[a2].x, vertices[a2].y, vertices[b1].x, vertices[b1].y, vertices[b2].x, vertices[b2].y );
                        // printf("%d\n", zz);
                        return "superposta";
                    }
                }
            }
        }
    }

    return NULL;
}

void print_vertices(Vertex *vertices, int n) {
    printf("VÉRTICES:\n");
    for (int i = 0; i < n; ++i) {
        printf("  v%d: (%d, %d)\n", i + 1, vertices[i].x, vertices[i].y);
    }
}

void print_faces(Face *faces, int n) {
    printf("\nFACES:\n");
    for (int i = 0; i < n; ++i) {
        printf("  f%d: ", faces[i].id + 1);
        HalfEdge *start = faces[i].outer_component;
        HalfEdge *e = start;
        do {
            printf("v%d ", e->origin->id + 1);
            e = e->next;
        } while (e != start);
        printf("\n");
    }
}


void printRawFaces(RawFace *RawFace, int n_faces)
{
    for(int i = 0 ; i < n_faces; i++)
    {
        for(int j = 0; j < RawFace[i].n_vertices; j++)
        {
            printf("%d ", RawFace[i].vertex_indices[j]);
        }

        printf("\n");
    }
}

void printCoords(Coord *coords, int n_vertices)
{
    for(int i = 0; i < n_vertices; i++)
    {
        printf("V%d: (%d, %d)\n", i+1, coords[i].x, coords[i].y );
    }
}

void set_incident_edges_by_angle(Vertex *vertices, int n_vertices, HalfEdge **all_edges, int edge_count) {
    for (int i = 0; i < n_vertices; ++i) {
        Vertex *v = &vertices[i];
        double min_angle = 3.14 * 2;  // Initialize with a large value
        HalfEdge *best = NULL;

        for (int j = 0; j < edge_count; ++j) {
            HalfEdge *he = all_edges[j];
            if (he->origin != v) continue;

            // Calculate vector from vertex to next vertex in half-edge
            int dx = he->next->origin->x - v->x;
            int dy = he->next->origin->y - v->y;

            // Calculate angle using dot product and magnitude
            double dot_product = v->x * dx + v->y * dy;
            double magnitude_v = sqrt(v->x * v->x + v->y * v->y);
            double magnitude_he = sqrt(dx * dx + dy * dy);
            double angle = acos(dot_product / (magnitude_v * magnitude_he));

            // Check if angle is smaller than current minimum
            if (angle < min_angle) {
                min_angle = angle;
                best = he;
            }
        }

        if (best) v->incident_edge = best;
    }
}

Face *build_dcel(Vertex *vertices, RawFace *raw_faces, int n_faces, HalfEdge **all_edges, int *edge_count, int n_vertices) {
    Face *faces = malloc(n_faces * sizeof(Face));
    int e_id = 1;  // IDs iniciam em 1

    for (int i = 0; i < n_faces; ++i) {
        faces[i].id = i + 1;
        faces[i].outer_component = NULL;
        int nv = raw_faces[i].n_vertices;
        int *idx = raw_faces[i].vertex_indices;
        
        // printf("face: %d, Vertices:%d\n", faces[i].id, nv);
        HalfEdge **face_edges = malloc(nv * sizeof(HalfEdge*));

        for (int j = 0; j < nv; ++j) {
            HalfEdge *he = malloc(sizeof(HalfEdge));
            he->id = e_id++;
            he->origin = &vertices[idx[j] - 1];
            he->face = &faces[i];
            he->twin = NULL;
            he->next = NULL;
            he->prev = NULL;
            face_edges[j] = he;
            // printf("Vertice id:%d  ",idx[j]);
            // printf("Aresta:%d Origem:%d Face:%d\n", he->id, he->origin->id, he->face->id);

            
            all_edges[*edge_count] = he;
            (*edge_count)++;
        }



        for (int j = 0; j < nv; ++j) {
            face_edges[j]->next = face_edges[(j + 1) % nv];
            face_edges[(j + 1) % nv]->prev = face_edges[j];
        }

        faces[i].outer_component = face_edges[0];
        free(face_edges);
    }

    // Conectar twins
    for (int i = 0; i < *edge_count; ++i) {
        HalfEdge *e1 = all_edges[i];
        if (e1->twin) continue;

        for (int j = i + 1; j < *edge_count; ++j) {
            HalfEdge *e2 = all_edges[j];
            if (e2->twin) continue;

            if (e1->origin == e2->next->origin &&
                e2->origin == e1->next->origin) {
                e1->twin = e2;
                e2->twin = e1;
                break;
            }
        }
    }

    int count = *edge_count;

    set_incident_edges_by_angle(vertices, n_vertices, all_edges, count);

    return faces;
}

void print_dcel_output(Vertex *vertices, int n_vertices,
                       Face *faces, int n_faces,
                       HalfEdge **edges, int edge_count) {
    int m = edge_count / 2;
    printf("%d %d %d\n", n_vertices, m, n_faces);

    // Imprime vértices
    for (int i = 0; i < n_vertices; ++i) {
        int he_id = vertices[i].incident_edge->id;
        printf("%d %d %d\n", vertices[i].x, vertices[i].y, he_id);
    }

    // Imprime faces
    for (int i = 0; i < n_faces; ++i) {
        printf("%d\n", faces[i].outer_component->id);
    }

    // Imprime semi-arestas
    for (int i = 0; i < edge_count; ++i) {
        HalfEdge *he = edges[i];
        int origin_id = he->origin->id + 1;
        int twin_id = he->twin->id;
        int face_id = he->face->id;
        int next_id = he->next->id;
        int prev_id = he->prev->id;
        printf("%d %d %d %d %d\n", origin_id, twin_id, face_id, next_id, prev_id);
    }
}


int main() {
    int n_vertices, n_faces;
    scanf("%d %d", &n_vertices, &n_faces);


    Coord *coords = read_vertices(n_vertices);
    RawFace *raw_faces = read_faces(n_faces);

    // Criar vértices
    Vertex *vertices = malloc(n_vertices * sizeof(Vertex));
    for (int i = 0; i < n_vertices; i++) {
        vertices[i].id = i;
        vertices[i].x = coords[i].x;
        vertices[i].y = coords[i].y;
        vertices[i].incident_edge = NULL;
    }

    const char *erro_topo = check_edge_counts(raw_faces, n_faces);
    if (erro_topo) {
        printf("%s\n", erro_topo);
        return 1;
    }

    // printCoords(coords, n_vertices);
    
    const char *erro_intersec = check_face_intersections(coords, raw_faces, n_faces);
    if (erro_intersec) {
        printf("%s\n", erro_intersec);
        return 1;
    }
    
    int max_half_edges = 2 * n_faces * 10;  // número máximo aproximado
    HalfEdge **all_edges = malloc(max_half_edges * sizeof(HalfEdge*));
    int edge_count = 0;
    
    // print_vertices(vertices, n_vertices);
    // printRawFaces(raw_faces, n_faces);

    Face *faces = build_dcel(vertices, raw_faces, n_faces, all_edges, &edge_count, n_vertices);
    
    print_dcel_output(vertices, n_vertices, faces, n_faces, all_edges, edge_count);

    for (int i = 0; i < n_faces; ++i)
        free(raw_faces[i].vertex_indices);
    free(coords);
    
    free(raw_faces);
    free(vertices);
    free(edge_map);
    free(edge_counts);
    free(faces);
    free(all_edges);
    return 0;
}

