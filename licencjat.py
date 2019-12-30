import IsoSpecPy


def main():

    class FlowGraph(object):
        """
          Dla podanej listy posortowanych krawędzi postaci [wierzchołek1,wierzchołek2,przepływ=0] utworzonej z widma
        eksperymentalnego tworzy graf, pozwala dodawać do niego nowe wierzchołki z widma teoretycznego za pomocą
        add_new_peak(). Dodajac wierzcholek znajduje najtanszy mozliwy przeplyw [ calculate_flow()].
        Funkcja pomocnicza flow_actualization() aktualizuje przepływy i ich koszt a edge_cost() wylicza koszt przeplywu
        na podstawie roznicy mas pików.
          """

        def __init__(self, *args, formula=None, edge_list=None):
            try:
                self.peakGenerator = iter(IsoSpecPy.IsoLayeredGenerator(formula=formula))
            except Exception as e:
                raise e
            self.dist = 0
            self.suma = 0.0
            self.edgeList = edge_list

        def edge_cost(self, index):
            return self.edgeList[index][1][0] - self.edgeList[index][0][0]

        def flow_actualization(self, start, stop):

            direction = 1
            if start >= stop:  # flow right to left
                start, stop, direction = stop, start, -1
            lef_end_flag = True if not start == 0 else False
            right_end_flag = True if not stop == len(self.edgeList) - 1 else False
            prob_difference = self.edgeList[stop][1][1] - self.edgeList[start][0][1]
            if prob_difference >= 0:
                flow = self.edgeList[start][0][1]
                self.edgeList[start][0] = (self.edgeList[start][0][0], 0)
                if lef_end_flag:
                    self.edgeList[start - 1][1] = (self.edgeList[start][0][0], 0)
                self.edgeList[stop][1] = (self.edgeList[stop][1][0], prob_difference)
                if right_end_flag:
                    self.edgeList[stop + 1][0] = (self.edgeList[stop][1][0], prob_difference)
            else:
                flow = self.edgeList[stop][1][1]
                self.edgeList[stop][1] = (self.edgeList[stop][1][0], 0)
                if right_end_flag:
                    self.edgeList[stop + 1][0] = (self.edgeList[stop][1][0], 0)
                self.edgeList[start][0] = (self.edgeList[start][0][0], -prob_difference)
                if lef_end_flag:
                    self.edgeList[start - 1][1] = (self.edgeList[start][0][0], -prob_difference)
            for x in range(start, stop + 1):
                self.dist -= abs(self.edgeList[x][2] * self.edge_cost(x))
                self.edgeList[x][2] += flow * direction
                self.dist += abs(self.edgeList[x][2] * self.edge_cost(x))

        def calculate_flow(self, edgeindex):
            flg = 1  # jeśli 0 to nowy wierzchołek jest skrajnie lewy
            if edgeindex == 0 or edgeindex == len(self.edgeList) - 1:
                left_index, right_index = edgeindex, edgeindex
            else:
                left_index, right_index = edgeindex, edgeindex + 1
            dist_sum_left = self.edge_cost(left_index)
            dist_sum_right = self.edge_cost(right_index)
            if right_index == 0:  # Jesli nowy pik jest skrajnie lewy nie mozna wydluzac lewej sciezki
                dist_sum_left = float('inf')
                flg = 0
            if left_index == len(self.edgeList) - 1:  # analogicznie jesli jest skrajnie prawy
                dist_sum_right = float('inf')
            while self.edgeList[edgeindex][flg][1] > 0 and \
                    (dist_sum_left != float('inf') or dist_sum_right != float('inf')):
                if dist_sum_left <= dist_sum_right:
                    if self.edgeList[left_index][0][1] > 0:
                        self.flow_actualization(start=edgeindex, stop=left_index)

                    else:
                        left_index -= 1
                        dist_sum_left += self.edge_cost(left_index) if left_index >= 0 else float('inf')
                else:
                    if self.edgeList[right_index][1][1] > 0:
                        self.flow_actualization(start=edgeindex + flg, stop=right_index)
                    else:
                        right_index += 1
                        dist_sum_right += self.edge_cost(right_index) if right_index <= len(
                            self.edgeList) - 1 else float('inf')
            if dist_sum_left == float('inf') and dist_sum_right == float('inf'):
                raise ValueError
                # narazie jakikolwiek error, moze rzucic w sytuacji
                # gdy widmo eksperymentalne sumuje sie do mniejszej wartosci niz generowane

        def add_new_peak(self):
            new_peak = next(self.peakGenerator)
            self.suma += new_peak[1]
            if self.edgeList[0][0][0] >= new_peak[0]:
                self.edgeList.insert(0, [new_peak, self.edgeList[0][0], 0])
                self.calculate_flow(0)
            elif self.edgeList[-1][0][0] <= new_peak[0]:
                self.edgeList.append([self.edgeList[-1][1], new_peak, 0])
                self.calculate_flow(len(self.edgeList) - 1)
            else:
                for counter, value in enumerate(self.edgeList):
                    if value[0][0] <= new_peak[0] <= value[1][0]:
                        new_edge1 = [value[0], new_peak, 0]
                        new_edge2 = [new_peak, value[1], 0]
                        self.edgeList.pop(counter)
                        self.edgeList.insert(counter, new_edge1)
                        self.edgeList.insert(counter + 1, new_edge2)
                        self.calculate_flow(counter)
                        break

# bardzo wstępne testy
    eksperymentalne = IsoSpecPy.IsoTotalProb(formula="C100H212N15O7", prob_to_cover=0.995)
    moja_lista = []
    for mass, prob in eksperymentalne:
        moja_lista.append((mass, prob))
    moja_lista = sorted(moja_lista, key=lambda tup: tup[0])
    lista2 = []
    for i in range(len(moja_lista) - 1):
        lista2.append([moja_lista[i], moja_lista[i + 1], 0])
    graf = FlowGraph(edge_list=lista2, formula="C101H213N16O7")
    while graf.suma < 0.99:
        graf.add_new_peak()
    print(graf.dist)


if __name__ == '__main__':
    main()
