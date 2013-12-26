CFLAGS=-Wall -O3

IndEvi: bin/IndEvi

bin/IndEvi: .obj/IndEvi.o .obj/PTNode.o .obj/TreeWriter.o .obj/BitDb.o .obj/Itemset.o .obj/Transactionset.o .obj/CartesianProduct.o .obj/CartesianProductDb.o .obj/IndEviStatic.o
	g++ $(CFLAGS) $^ -o bin/IndEvi

.obj/%.o: src/%.cpp
	g++ $(CFLAGS) $< -c -o $@

IndEviRe: bin/IndEviRe	

bin/IndEviRe: .obj/IndEviRe.o .obj/PTNode.o .obj/TreeWriter.o .obj/BitDb.o .obj/Itemset.o .obj/Transactionset.o .obj/CartesianProduct.o .obj/CartesianProductDb.o .obj/IndEviStatic.o
	g++ $(CFLAGS) .obj/PTNode.o .obj/TreeWriter.o .obj/BitDb.o .obj/Itemset.o .obj/Transactionset.o .obj/CartesianProduct.o .obj/CartesianProductDb.o .obj/IndEviStatic.o .obj/IndEviRe.o -o bin/IndEviRe	
	
	
clean:
	rm -f .obj/*.o bin/IndEvi bin/postprocessing bin/IndEviRe
