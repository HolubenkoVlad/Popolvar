#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/timeb.h>

#define PARENT(x) (x-1)/2
#define LEFT(x) 2*x+1
#define RIGHT(x) 2*i+2

typedef struct HeapNode {
	int x;
	int y;
	int weight;
	int prev_x;
	int prev_y;
}HeapNode;

typedef struct MinHeap {
	HeapNode** data;
	int heap_size;
}MinHeap;

typedef struct Node {
	int weight;
	int visited;
	int prev_coords[2];
}Node;


int* zachran_princezne(char**, int, int, int, int*);

int* dijkstra(Node**, int, int, int, int, int, int, int*, int*);


Node* createNode(int weight) {
	Node* node = (Node*)malloc(sizeof(Node));
	if (node == NULL) {
		printf("can't alloc the memory\n");
		return NULL;
	}
	node->weight = weight;
	node->visited = 0;
	return node;
}

HeapNode* createHeapNode(int x, int y, int weight, int i, int j) {
	HeapNode* heapNode = (HeapNode*)malloc(sizeof(HeapNode));
	heapNode->x = x;
	heapNode->y = y;
	heapNode->prev_x = i;
	heapNode->prev_y = j;
	heapNode->weight = weight;
	return heapNode;
}

MinHeap* createMinHeap(int size) {
	MinHeap* min_heap = (MinHeap*)malloc(sizeof(MinHeap));
	min_heap->data = (HeapNode**)calloc(size, sizeof(HeapNode*));
	min_heap->heap_size = 0;
	return min_heap;
}

void swap(HeapNode** x, HeapNode** y)
{
	HeapNode* temp = *x;
	*x = *y;
	*y = temp;
}

void insertNode(HeapNode* heapNode, MinHeap* heap)
{
	heap->heap_size++;
	int i = heap->heap_size - 1;
	heap->data[i] = heapNode;
	while (i != 0 && heap->data[PARENT(i)]->weight > heap->data[i]->weight)
	{
		swap(&heap->data[i], &heap->data[PARENT(i)]);
		i = PARENT(i);
	}
}


void MinHeapify(int i, MinHeap* heap)
{
	int l = LEFT(i);
	int r = RIGHT(i);
	int smallest = i;
	if (l < heap->heap_size && heap->data[l]->weight < heap->data[i]->weight)
		smallest = l;
	if (r < heap->heap_size && heap->data[r]->weight < heap->data[smallest]->weight)
		smallest = r;
	if (smallest != i)
	{
		swap(&heap->data[i], &heap->data[smallest]);
		MinHeapify(smallest, heap);
	}
}

HeapNode* getMin(MinHeap* heap)
{
	if (heap->heap_size <= 0)
		return NULL;
	if (heap->heap_size == 1)
	{
		heap->heap_size--;
		return heap->data[0];
	}
	HeapNode* root = heap->data[0];
	heap->data[0] = heap->data[heap->heap_size - 1];
	heap->heap_size--;
	MinHeapify(0, heap);
	return root;
}

float calculate_time(struct timeb start, struct timeb end) {
	return (float)(((end.time - start.time) * 1000) + end.millitm - start.millitm) / 1000;
}


int main() {
	char** mapa = NULL;
	int i, test, dlzka_cesty, cas = 0, * cesta = NULL;
	struct timeb start, end;
	int n = 0, m = 0, t = 0;
	FILE* f;
	while (1) {
		printf("Zadajte cislo testu (0 ukonci program):\n");
		scanf("%d", &test);
		dlzka_cesty = 0;
		n = m = t = 0;
		switch (test) {
		case 0://ukonci program
			return 0;
		case 1://nacitanie mapy zo suboru
			f = fopen("test.txt", "r");
			if (f)
				fscanf(f, "%d %d %d", &n, &m, &t);
			else
				return;
			mapa = (char**)malloc(n * sizeof(char*));
			for (i = 0; i < n; i++) {
				mapa[i] = (char*)malloc(m * sizeof(char));
				for (int j = 0; j < m; j++) {
					char policko = fgetc(f);
					if (policko == '\n') policko = fgetc(f);
					mapa[i][j] = policko;
				}
			}
			fclose(f);
			ftime(&start);
			cesta = zachran_princezne(mapa, n, m, t, &dlzka_cesty);
			ftime(&end);
			if(cesta!=NULL)
			printf("cas vykonavania: %.3f s\n", calculate_time(start, end));
			break;
		case 2://nacitanie preddefinovanej mapy
			n = 10;
			m = 10;
			t = 12;
			mapa = (char**)malloc(n * sizeof(char*));
			mapa[0] = "CCHCNHCCHN";
			mapa[1] = "NNCCCHHCCC";
			mapa[2] = "DNCCNNHHHC";
			mapa[3] = "CHHHCCCCCC";
			mapa[4] = "CCCCCNHHHH";
			mapa[5] = "PCHCCCNNNN";
			mapa[6] = "NNNNNHCCCC";
			mapa[7] = "CCCCCPCCCC";
			mapa[8] = "CCCNNHHHHH";
			mapa[9] = "HHHPCCCCCC";
			ftime(&start);
			cesta = zachran_princezne(mapa, n, m, t, &dlzka_cesty);
			ftime(&end);
			if (cesta != NULL)
			printf("cas vykonavania: %.3f s\n", calculate_time(start, end));
			break;
		}
		if (cesta == NULL) {
			return 0;
		}
		for (int i = 0; i < dlzka_cesty; i++) {
			printf("%d %d\n", cesta[i * 2 + 1], cesta[i * 2]);
			if (mapa[cesta[i * 2]][cesta[i * 2 + 1]] == 'H')
				cas += 2;
			else
				cas += 1;
			if (mapa[cesta[i * 2]][cesta[i * 2 + 1]] == 'D' && cas > t)
				printf("Nestihol si zabit draka!\n");
			if (mapa[cesta[i * 2]][cesta[i * 2 + 1]] == 'N')
				printf("Prechod cez nepriechodnu prekazku!\n");
			if (i > 0 && abs(cesta[i * 2 + 1] - cesta[(i - 1) * 2 + 1]) + abs(cesta[i * 2] - cesta[(i - 1) * 2]) > 1)
				printf("Neplatny posun Popolvara!\n");

		}
		printf("cas cesty: %d\n", cas);
		free(cesta);
		cas = 0;
		if (test == 1)
			for (int i = 0; i < n; i++) {
				free(mapa[i]);
			}
		free(mapa);
	}
	return 0;
}

int map_value(char symbol) {
	if (symbol == 'C' || symbol == 'c' || symbol == 'D' || symbol == 'P')
		return 1;
	else if (symbol == 'H' || symbol == 'h')
		return 2;
	else return -1;
}

//vytvorim mapu s uloženou veľkosťou a súradnicami predchdzajucej bodky
Node** createDependencyMap(char** map, int n, int m, int** persons) {
	*persons = (int*)calloc((6 * 2) + 1, sizeof(int));
	Node** dependency_map = (Node**)malloc(n * m * sizeof(Node*));
	int count_princ = 1, dragon_counter = 0;
	int flag1 = 1, flag2 = 1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			*(dependency_map + i * m + j) = createNode(map_value(map[i][j]));
			if (map[i][j] == 'P' || map[i][j] == 'p') {
				(*persons)[count_princ * 2 + 1] = i;
				(*persons)[count_princ * 2 + 2] = j;
				count_princ++;
				flag1 = 0;
			}
			else if (map[i][j] == 'D' || map[i][j] == 'd') {
				(*persons)[1] = i;
				(*persons)[2] = j;
				flag2 = 0;
				dragon_counter++;
			}
		}
	}
	if (count_princ > 6 || dragon_counter > 1)
	{
		printf("Viac ako 5 princeznej alebo viac ako 1 drak!\n");
		return NULL;
	}
	(*persons)[0] = count_princ;
	if (flag1 == 1) {
		printf("Neexistuje ziadna princezna!\n");
		return NULL;
	}
	else if (flag2 == 1) {
		printf("Neexistuje ziaden drak!\n");
		return NULL;
	}
	else return dependency_map;
}

int* dijkstra(Node** dependency_map, int n, int m, int startX, int startY, int endX, int endY, int* dlzka_cesty, int* way_weight) {
	int i = startX, j = startY;
	MinHeap* minHeap = createMinHeap(n * m);
	HeapNode* temp;
	while (1) {
		if (i == endX && j == endY)
			break;
		for (int k = j - 1, p = i - 1; k <= j + 1; k += 2, p += 2) {
			//Pozerám sa na unvisited body vľavo a vpravo od bodu a pridavam do min heap
			if ((k >= 0 && k < m) && ((*(dependency_map + i * m + k))->weight >= 0 && (*(dependency_map + i * m + k))->visited == 0)) {
				(*(dependency_map + i * m + k))->visited = 1;
				insertNode(createHeapNode(i, k, (*(dependency_map + i * m + j))->weight + (*(dependency_map + i * m + k))->weight, i, j), minHeap);
			}
			//Pozerám sa na unvisited body hore a dole od bodu a pridavam do min heap
			if ((p >= 0 && p < n) && ((*(dependency_map + p * m + j))->weight >= 0 && (*(dependency_map + p * m + j))->visited == 0)) {
				(*(dependency_map + p * m + j))->visited = 1;
				insertNode(createHeapNode(p, j, (*(dependency_map + i * m + j))->weight + (*(dependency_map + p * m + j))->weight, i, j), minHeap);
			}
		}
		//vyberam minimalny bod
		temp = getMin(minHeap);
		if (temp == NULL)
			return NULL;
		//aktualizujem veľkosť bodky a koordinaty predchadzajucej bodky
		(*(dependency_map + temp->x * m + temp->y))->weight = temp->weight;
		(*(dependency_map + temp->x * m + temp->y))->prev_coords[0] = temp->prev_x;
		(*(dependency_map + temp->x * m + temp->y))->prev_coords[1] = temp->prev_y;
		i = temp->x;
		j = temp->y;
	}
	//najdem dlžku cesty
	int temp1, temp2, counter = 0;
	for (int i = endX, j = endY; ; i = temp1, j = temp2) {
		temp1 = (*(dependency_map + i * m + j))->prev_coords[0];
		temp2 = (*(dependency_map + i * m + j))->prev_coords[1];
		counter++;
		//printf("%d %d\n", temp1, temp2);
		if (i == startX && j == startY)
			break;
	}
	//for (int i = 0; i < minHeap->heap_size; i++)
	//	free(minHeap->data[i]);
	free(minHeap);
	//printf("\n");
	counter *= 2;
	*dlzka_cesty = counter;
	*way_weight = (*(dependency_map + endX * m + endY))->weight;
	//vytvorim cestu
	int* cesta = (int*)calloc(counter, sizeof(int));
	cesta[counter - 2] = endX; cesta[counter - 1] = endY;
	counter -= 3;
	for (int i = endX, j = endY; counter > 0; i = temp1, j = temp2, counter -= 2) {
		temp1 = (*(dependency_map + i * m + j))->prev_coords[0];
		temp2 = (*(dependency_map + i * m + j))->prev_coords[1];
		cesta[counter - 1] = temp1;
		cesta[counter] = temp2;
	}
	return cesta;
}

int factorial(int number)
{
	if (number == 1)return number;
	return factorial(number - 1) * number;
}

void permut_swap(int* a, int i, int j)
{
	int s = a[i];
	a[i] = a[j];
	a[j] = s;
}
//nachadza všetci permutacie
int PermutationSet(int* perm_array, int number)
{
	int j = number - 2;
	while (j != -1 && perm_array[j] >= perm_array[j + 1]) j--;
	if (j == -1)
		return 0;
	int k = number - 1;
	while (perm_array[j] >= perm_array[k]) k--;
	permut_swap(perm_array, j, k);
	int l = j + 1, r = number - 1;
	while (l < r)
		permut_swap(perm_array, l++, r--);
	return 1;
}

void free_map(Node** node, int n, int m) {
	for (int i = 0; i < n * m; i++) {
		free(*(node + i));
	}
	free(node);
}
//vytvara copiu matice s uloženou veľkosťou a súradnicami predchdzajucej bodky
Node** full_copy(Node** first, Node** second, int n, int m) {
	//if (second != NULL)
	//{
		//free_map(second, n, m);
		//free(second);
	//}
	Node** new_second = (Node**)malloc(n * m * sizeof(Node*));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			*(new_second + i * m + j) = createNode((*(first + i * m + j))->weight);
		}
	}
	return new_second;
}

int* zachran_princezne(char** mapa, int n, int m, int t, int* dlzka_cesty) {
	int* persons_coord;
	int way_weight[1];
	//vytvorim maticu s velkostou a predchadzajucej bodkou, a najdem všetkych hrdinov
	Node** dependency_map = createDependencyMap(mapa, n, m, &persons_coord);
	if (dependency_map == NULL) {
		return NULL;
	}
	//------------------
	//vytvorim copiu tej matice, aby pracovat s nou v dijkstra
	Node** temp_dep_map = NULL;
	temp_dep_map = full_copy(dependency_map, temp_dep_map, n, m);
	//------------------
	//nachadzam cestu
	int* cesta = dijkstra(temp_dep_map, n, m, 0, 0, persons_coord[1], persons_coord[2], dlzka_cesty, way_weight);
	if (cesta == NULL) {
		printf("Nemozem dosiahnut hrdinu!\n");
		return NULL;
	}
	//------------------
	//najdem všetci permutacie
	int fact = factorial(persons_coord[0] - 1);
	int** permutation_array = (int**)malloc(fact * sizeof(int*));
	int count_of_princess = persons_coord[0] - 1;
	for (int i = 0; i < fact; i++)
		permutation_array[i] = (int*)calloc(count_of_princess, sizeof(int));
	int* a = (int*)calloc(count_of_princess, sizeof(int));
	for (int i = 0; i < count_of_princess; i++)
		a[i] = i + 1;
	int flag, counter = 0;
	do {
		memcpy(permutation_array[counter], a, ((persons_coord[0] - 1) * 2) * sizeof(int));
		flag = PermutationSet(a, persons_coord[0] - 1);
		counter++;
	} while (flag != 0);
	//---------------------
	//hľadanie najlepšej cesty
	int* prev_cesta = NULL, length[1], length_prev[1], weight[1], weight_prev[1], * temp_cesta = NULL, length_temp[1], weight_temp[1],
		* first_cesta = (int*)calloc(dlzka_cesty[0], sizeof(int));
	int j = 0;
	*length_temp = 0;
	//uchovam cestu od hrdinu do draka
	memcpy(first_cesta, cesta, (dlzka_cesty[0]) * sizeof(int));
	//---------------
	for (int i = 0; i < fact; i++) {
		if (prev_cesta != NULL)
			free(prev_cesta);
		temp_dep_map = full_copy(dependency_map, temp_dep_map, n, m);
		//najsť cestu od draka do prvej princezne v zozname
		prev_cesta = dijkstra(temp_dep_map, n, m, persons_coord[1], persons_coord[2], persons_coord[permutation_array[i][j] * 2 + 1], persons_coord[permutation_array[i][j] * 2 + 2], length, weight);
		if (prev_cesta == NULL) {
			printf("Nemozem dosiahnut hrdinu\n");
			return NULL;
		}
		//spojiť cesty
		cesta = (int*)realloc(cesta, (dlzka_cesty[0] + (length[0]) - 2) * sizeof(int));
		memmove(cesta + dlzka_cesty[0] - 2, prev_cesta, (length[0]) * sizeof(int));
		weight[0] = way_weight[0] + weight[0] - 1;
		length[0] = (dlzka_cesty[0] + length[0] - 2);
		//-----------
		for (j = 1; j < persons_coord[0] - 1; j++) {
			length_prev[0] = length[0];
			weight_prev[0] = weight[0];
			if (prev_cesta != NULL)
				free(prev_cesta);
			temp_dep_map = full_copy(dependency_map, temp_dep_map, n, m);
			//najsť cestu od prvej princezne v zozname do dalšej
			prev_cesta = dijkstra(temp_dep_map, n, m, persons_coord[permutation_array[i][j - 1] * 2 + 1], persons_coord[permutation_array[i][j - 1] * 2 + 2], persons_coord[permutation_array[i][j] * 2 + 1], persons_coord[permutation_array[i][j] * 2 + 2], length, weight);
			if (prev_cesta == NULL) {
				printf("Nemozem dosiahnut hrdinu\n");
				return NULL;
			}
			//spojit cesty
			cesta = (int*)realloc(cesta, (length_prev[0] + (length[0]) - 2) * sizeof(int));
			memmove(cesta + length_prev[0] - 2, prev_cesta, (length[0]) * sizeof(int));
			weight[0] = weight[0] + weight_prev[0] - 1;
			length[0] = (length_prev[0] + length[0] - 2);
			//-------------
		}
		//ak cesta ma mensiu veľkost než predchádzajúca cesta uchovať cestu
		if (temp_cesta == NULL || weight[0] < weight_temp[0]) {
			if (temp_cesta != NULL)
				free(temp_cesta);
			temp_cesta = (int*)calloc(length[0], sizeof(int));
			memmove(temp_cesta, cesta, length[0] * sizeof(int));
			length_temp[0] = length[0];
			weight_temp[0] = weight[0];
		}
		free(cesta);
		cesta = (int*)calloc(*dlzka_cesty, sizeof(int));
		//obnoviť cestu od hrdinu k draku
		memmove(cesta, first_cesta, (dlzka_cesty[0]) * sizeof(int));
		j = 0;
	}
	free_map(dependency_map, n, m);
	free_map(temp_dep_map, n, m);
	free(first_cesta);
	free(prev_cesta);
	free(cesta);
	//free(a);
	for (int i = 0; i < fact; i++) {
		//free(permutation_array[i]);
	}
	free(permutation_array);
	*dlzka_cesty = length_temp[0] / 2;
	return temp_cesta;
}