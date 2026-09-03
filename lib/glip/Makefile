
LIB=./lib
INCLUDE=./include
BIN=./bin
SRC=./src
OBJ=./obj
DEMO=./demo

#Compiladores
CXX=g++

#-march=native
#FLAGS= -Wall -msse -fopenmp
#FLAGS= -Wall -O3 -msse -fopenmp
#FLAGS= -Wall -g
#FLAGS= -Wall -g -O3
FLAGS= -Wall -O3

LINKS= -lpthread -lz -lm 


all: libglip $(BIN)/areaclosing $(BIN)/volumeclosing $(BIN)/superpixels $(BIN)/supervoxels
	@echo
	@echo "libglip.a built..."
	@echo "DONE."

libglip: $(LIB)/libglip.a

GLIPLIB  = -L$(LIB) -lglip
GLIPFLAGS  = -I$(INCLUDE)


$(LIB)/libglip.a: \
$(OBJ)/glip_common.o \
$(OBJ)/glip_adjrel.o \
$(OBJ)/glip_image32.o \
$(OBJ)/glip_image32f.o \
$(OBJ)/glip_image64.o \
$(OBJ)/glip_image8.o \
$(OBJ)/glip_cimage.o \
$(OBJ)/glip_cimage32f.o \
$(OBJ)/glip_morphology.o \
$(OBJ)/glip_morphology3.o \
$(OBJ)/glip_queue.o \
$(OBJ)/glip_gpqueue_by_Falcao.o \
$(OBJ)/glip_graph.o \
$(OBJ)/glip_imagegraph.o \
$(OBJ)/glip_imagegraphpx.o \
$(OBJ)/glip_layeredgraph.o \
$(OBJ)/glip_heap.o \
$(OBJ)/glip_heap32.o \
$(OBJ)/glip_heap64f.o \
$(OBJ)/glip_heap_lex.o \
$(OBJ)/glip_heap32fi_lex.o \
$(OBJ)/glip_heappair.o \
$(OBJ)/glip_ift.o \
$(OBJ)/glip_string.o \
$(OBJ)/glip_arraylist.o \
$(OBJ)/glip_filelist.o \
$(OBJ)/glip_color.o \
$(OBJ)/glip_highlight.o \
$(OBJ)/glip_superpixels.o \
$(OBJ)/glip_evaluation.o \
$(OBJ)/glip_scene.o \
$(OBJ)/glip_scene32.o \
$(OBJ)/glip_scene32f.o \
$(OBJ)/glip_scene64.o \
$(OBJ)/glip_scene64f.o \
$(OBJ)/glip_scene16.o \
$(OBJ)/glip_scene8.o \
$(OBJ)/glip_adjrel3.o \
$(OBJ)/glip_scnmath.o \
$(OBJ)/glip_matrix.o \
$(OBJ)/glip_bmap.o \
$(OBJ)/glip_set.o \
$(OBJ)/glip_stack.o \
$(OBJ)/glip_analysis.o \
$(OBJ)/glip_analysis3.o \
$(OBJ)/glip_pqueue32.o \
$(OBJ)/glip_pqueue16.o \
$(OBJ)/glip_pqueue8.o \
$(OBJ)/glip_hashtable.o \
$(OBJ)/glip_curve.o \
$(OBJ)/glip_radiometric.o \
$(OBJ)/glip_radiometric3.o \
$(OBJ)/glip_filtering.o \
$(OBJ)/glip_filtering3.o \
$(OBJ)/glip_segmentation3.o \
$(OBJ)/glip_gradient3.o \
$(OBJ)/glip_marker.o \
$(OBJ)/glip_marker3.o \
$(OBJ)/glip_features.o \
$(OBJ)/glip_graphclustering.o \
$(OBJ)/glip_attributelist.o \
$(OBJ)/glip_configfile.o \
$(OBJ)/glip_geometry.o \
$(OBJ)/glip_ctree.o \
$(OBJ)/glip_llw.o \
$(OBJ)/glip_draw.o \
$(OBJ)/glip_band.o \
$(OBJ)/glip_dicom.o \
$(OBJ)/glip_bzlib.o \
$(OBJ)/nifti1_io.o \
$(OBJ)/znzlib.o
	ar csr $(LIB)/libglip.a \
$(OBJ)/glip_common.o \
$(OBJ)/glip_adjrel.o \
$(OBJ)/glip_image32.o \
$(OBJ)/glip_image32f.o \
$(OBJ)/glip_image64.o \
$(OBJ)/glip_image8.o \
$(OBJ)/glip_cimage.o \
$(OBJ)/glip_cimage32f.o \
$(OBJ)/glip_morphology.o \
$(OBJ)/glip_morphology3.o \
$(OBJ)/glip_queue.o \
$(OBJ)/glip_gpqueue_by_Falcao.o \
$(OBJ)/glip_graph.o \
$(OBJ)/glip_imagegraph.o \
$(OBJ)/glip_imagegraphpx.o \
$(OBJ)/glip_layeredgraph.o \
$(OBJ)/glip_heap.o \
$(OBJ)/glip_heap32.o \
$(OBJ)/glip_heap64f.o \
$(OBJ)/glip_heap_lex.o \
$(OBJ)/glip_heap32fi_lex.o \
$(OBJ)/glip_heappair.o \
$(OBJ)/glip_ift.o \
$(OBJ)/glip_string.o \
$(OBJ)/glip_arraylist.o \
$(OBJ)/glip_filelist.o \
$(OBJ)/glip_color.o \
$(OBJ)/glip_highlight.o \
$(OBJ)/glip_superpixels.o \
$(OBJ)/glip_evaluation.o \
$(OBJ)/glip_scene.o \
$(OBJ)/glip_scene32.o \
$(OBJ)/glip_scene32f.o \
$(OBJ)/glip_scene64.o \
$(OBJ)/glip_scene64f.o \
$(OBJ)/glip_scene16.o \
$(OBJ)/glip_scene8.o \
$(OBJ)/glip_adjrel3.o \
$(OBJ)/glip_scnmath.o \
$(OBJ)/glip_matrix.o \
$(OBJ)/glip_bmap.o \
$(OBJ)/glip_set.o \
$(OBJ)/glip_stack.o \
$(OBJ)/glip_analysis.o \
$(OBJ)/glip_analysis3.o \
$(OBJ)/glip_pqueue32.o \
$(OBJ)/glip_pqueue16.o \
$(OBJ)/glip_pqueue8.o \
$(OBJ)/glip_hashtable.o \
$(OBJ)/glip_curve.o \
$(OBJ)/glip_radiometric.o \
$(OBJ)/glip_radiometric3.o \
$(OBJ)/glip_filtering.o \
$(OBJ)/glip_filtering3.o \
$(OBJ)/glip_segmentation3.o \
$(OBJ)/glip_gradient3.o \
$(OBJ)/glip_marker.o \
$(OBJ)/glip_marker3.o \
$(OBJ)/glip_features.o \
$(OBJ)/glip_graphclustering.o \
$(OBJ)/glip_attributelist.o \
$(OBJ)/glip_configfile.o \
$(OBJ)/glip_geometry.o \
$(OBJ)/glip_ctree.o \
$(OBJ)/glip_llw.o \
$(OBJ)/glip_draw.o \
$(OBJ)/glip_band.o \
$(OBJ)/glip_dicom.o \
$(OBJ)/glip_bzlib.o \
$(OBJ)/nifti1_io.o \
$(OBJ)/znzlib.o


$(OBJ)/glip_common.o: $(SRC)/glip_common.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_common.cpp -I$(INCLUDE) -o $(OBJ)/glip_common.o

$(OBJ)/glip_adjrel.o: $(SRC)/glip_adjrel.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_adjrel.cpp -I$(INCLUDE) -o $(OBJ)/glip_adjrel.o

$(OBJ)/glip_image32.o: $(SRC)/glip_image32.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_image32.cpp -I$(INCLUDE) -o $(OBJ)/glip_image32.o

$(OBJ)/glip_image32f.o: $(SRC)/glip_image32f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_image32f.cpp -I$(INCLUDE) -o $(OBJ)/glip_image32f.o

$(OBJ)/glip_image64.o: $(SRC)/glip_image64.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_image64.cpp -I$(INCLUDE) -o $(OBJ)/glip_image64.o

$(OBJ)/glip_image8.o: $(SRC)/glip_image8.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_image8.cpp -I$(INCLUDE) -o $(OBJ)/glip_image8.o

$(OBJ)/glip_cimage.o: $(SRC)/glip_cimage.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_cimage.cpp -I$(INCLUDE) -o $(OBJ)/glip_cimage.o

$(OBJ)/glip_cimage32f.o: $(SRC)/glip_cimage32f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_cimage32f.cpp -I$(INCLUDE) -o $(OBJ)/glip_cimage32f.o

$(OBJ)/glip_morphology.o: $(SRC)/glip_morphology.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_morphology.cpp -I$(INCLUDE) -o $(OBJ)/glip_morphology.o

$(OBJ)/glip_morphology3.o: $(SRC)/glip_morphology3.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_morphology3.cpp -I$(INCLUDE) -o $(OBJ)/glip_morphology3.o

$(OBJ)/glip_queue.o: $(SRC)/glip_queue.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_queue.cpp -I$(INCLUDE) -o $(OBJ)/glip_queue.o

$(OBJ)/glip_gpqueue_by_Falcao.o: $(SRC)/glip_gpqueue_by_Falcao.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_gpqueue_by_Falcao.cpp -I$(INCLUDE) -o $(OBJ)/glip_gpqueue_by_Falcao.o

$(OBJ)/glip_graph.o: $(SRC)/glip_graph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_graph.cpp -I$(INCLUDE) -o $(OBJ)/glip_graph.o

$(OBJ)/glip_graphclustering.o: $(SRC)/glip_graphclustering.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_graphclustering.cpp -I$(INCLUDE) -o $(OBJ)/glip_graphclustering.o

$(OBJ)/glip_imagegraph.o: $(SRC)/glip_imagegraph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_imagegraph.cpp -I$(INCLUDE) -o $(OBJ)/glip_imagegraph.o

$(OBJ)/glip_imagegraphpx.o: $(SRC)/glip_imagegraphpx.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_imagegraphpx.cpp -I$(INCLUDE) -o $(OBJ)/glip_imagegraphpx.o

$(OBJ)/glip_layeredgraph.o: $(SRC)/glip_layeredgraph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_layeredgraph.cpp -I$(INCLUDE) -o $(OBJ)/glip_layeredgraph.o

$(OBJ)/glip_heap.o: $(SRC)/glip_heap.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_heap.cpp -I$(INCLUDE) -o $(OBJ)/glip_heap.o

$(OBJ)/glip_heap32.o: $(SRC)/glip_heap32.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_heap32.cpp -I$(INCLUDE) -o $(OBJ)/glip_heap32.o

$(OBJ)/glip_heap64f.o: $(SRC)/glip_heap64f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_heap64f.cpp -I$(INCLUDE) -o $(OBJ)/glip_heap64f.o

$(OBJ)/glip_heap_lex.o: $(SRC)/glip_heap_lex.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_heap_lex.cpp -I$(INCLUDE) -o $(OBJ)/glip_heap_lex.o

$(OBJ)/glip_heap32fi_lex.o: $(SRC)/glip_heap32fi_lex.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_heap32fi_lex.cpp -I$(INCLUDE) -o $(OBJ)/glip_heap32fi_lex.o

$(OBJ)/glip_heappair.o: $(SRC)/glip_heappair.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_heappair.cpp -I$(INCLUDE) -o $(OBJ)/glip_heappair.o

$(OBJ)/glip_ift.o: $(SRC)/glip_ift.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_ift.cpp -I$(INCLUDE) -o $(OBJ)/glip_ift.o

$(OBJ)/glip_string.o: $(SRC)/glip_string.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_string.cpp -I$(INCLUDE) -o $(OBJ)/glip_string.o

$(OBJ)/glip_arraylist.o: $(SRC)/glip_arraylist.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_arraylist.cpp -I$(INCLUDE) -o $(OBJ)/glip_arraylist.o

$(OBJ)/glip_filelist.o: $(SRC)/glip_filelist.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_filelist.cpp -I$(INCLUDE) -o $(OBJ)/glip_filelist.o

$(OBJ)/glip_highlight.o: $(SRC)/glip_highlight.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_highlight.cpp -I$(INCLUDE) -o $(OBJ)/glip_highlight.o

$(OBJ)/glip_color.o: $(SRC)/glip_color.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_color.cpp -I$(INCLUDE) -o $(OBJ)/glip_color.o


$(OBJ)/glip_superpixels.o: $(SRC)/glip_superpixels.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_superpixels.cpp -I$(INCLUDE) -o $(OBJ)/glip_superpixels.o


$(OBJ)/glip_evaluation.o: $(SRC)/glip_evaluation.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_evaluation.cpp -I$(INCLUDE) -o $(OBJ)/glip_evaluation.o

$(OBJ)/glip_marker.o: $(SRC)/glip_marker.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_marker.cpp -I$(INCLUDE) -o $(OBJ)/glip_marker.o

$(OBJ)/glip_marker3.o: $(SRC)/glip_marker3.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_marker3.cpp -I$(INCLUDE) -o $(OBJ)/glip_marker3.o

$(OBJ)/glip_scene.o: $(SRC)/glip_scene.cpp $(INCLUDE)/glip_scene.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_scene.cpp -I$(INCLUDE) -o $(OBJ)/glip_scene.o

$(OBJ)/glip_scene32.o: $(SRC)/glip_scene32.cpp $(INCLUDE)/glip_scene32.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_scene32.cpp -I$(INCLUDE) -o $(OBJ)/glip_scene32.o

$(OBJ)/glip_scene64.o: $(SRC)/glip_scene64.cpp $(INCLUDE)/glip_scene64.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_scene64.cpp -I$(INCLUDE) -o $(OBJ)/glip_scene64.o

$(OBJ)/glip_scene32f.o: $(SRC)/glip_scene32f.cpp $(INCLUDE)/glip_scene32f.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_scene32f.cpp -I$(INCLUDE) -o $(OBJ)/glip_scene32f.o

$(OBJ)/glip_scene64f.o: $(SRC)/glip_scene64f.cpp $(INCLUDE)/glip_scene64f.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_scene64f.cpp -I$(INCLUDE) -o $(OBJ)/glip_scene64f.o

$(OBJ)/glip_scene16.o: $(SRC)/glip_scene16.cpp $(INCLUDE)/glip_scene16.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_scene16.cpp -I$(INCLUDE) -o $(OBJ)/glip_scene16.o

$(OBJ)/glip_scene8.o: $(SRC)/glip_scene8.cpp $(INCLUDE)/glip_scene8.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_scene8.cpp -I$(INCLUDE) -o $(OBJ)/glip_scene8.o

$(OBJ)/glip_adjrel3.o: $(SRC)/glip_adjrel3.cpp $(INCLUDE)/glip_adjrel3.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_adjrel3.cpp -I$(INCLUDE) -o $(OBJ)/glip_adjrel3.o

$(OBJ)/glip_scnmath.o: $(SRC)/glip_scnmath.cpp $(INCLUDE)/glip_scnmath.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_scnmath.cpp -I$(INCLUDE) -o $(OBJ)/glip_scnmath.o

$(OBJ)/glip_matrix.o: $(SRC)/glip_matrix.cpp $(INCLUDE)/glip_matrix.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_matrix.cpp -I$(INCLUDE) -o $(OBJ)/glip_matrix.o

$(OBJ)/glip_bmap.o: $(SRC)/glip_bmap.cpp $(INCLUDE)/glip_bmap.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_bmap.cpp -I$(INCLUDE) -o $(OBJ)/glip_bmap.o

$(OBJ)/glip_set.o: $(SRC)/glip_set.cpp $(INCLUDE)/glip_set.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_set.cpp -I$(INCLUDE) -o $(OBJ)/glip_set.o

$(OBJ)/glip_stack.o: $(SRC)/glip_stack.cpp $(INCLUDE)/glip_stack.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_stack.cpp -I$(INCLUDE) -o $(OBJ)/glip_stack.o

$(OBJ)/glip_analysis.o: $(SRC)/glip_analysis.cpp $(INCLUDE)/glip_analysis.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_analysis.cpp -I$(INCLUDE) -o $(OBJ)/glip_analysis.o

$(OBJ)/glip_analysis3.o: $(SRC)/glip_analysis3.cpp $(INCLUDE)/glip_analysis3.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_analysis3.cpp -I$(INCLUDE) -o $(OBJ)/glip_analysis3.o

$(OBJ)/glip_pqueue32.o: $(SRC)/glip_pqueue32.cpp $(INCLUDE)/glip_pqueue32.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_pqueue32.cpp -I$(INCLUDE) -o $(OBJ)/glip_pqueue32.o

$(OBJ)/glip_pqueue16.o: $(SRC)/glip_pqueue16.cpp $(INCLUDE)/glip_pqueue16.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_pqueue16.cpp -I$(INCLUDE) -o $(OBJ)/glip_pqueue16.o

$(OBJ)/glip_pqueue8.o: $(SRC)/glip_pqueue8.cpp $(INCLUDE)/glip_pqueue8.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_pqueue8.cpp -I$(INCLUDE) -o $(OBJ)/glip_pqueue8.o

$(OBJ)/glip_hashtable.o: $(SRC)/glip_hashtable.cpp $(INCLUDE)/glip_hashtable.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_hashtable.cpp -I$(INCLUDE) -o $(OBJ)/glip_hashtable.o

$(OBJ)/glip_curve.o: $(SRC)/glip_curve.cpp $(INCLUDE)/glip_curve.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_curve.cpp -I$(INCLUDE) -o $(OBJ)/glip_curve.o

$(OBJ)/glip_radiometric.o: $(SRC)/glip_radiometric.cpp $(INCLUDE)/glip_radiometric.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_radiometric.cpp -I$(INCLUDE) -o $(OBJ)/glip_radiometric.o

$(OBJ)/glip_radiometric3.o: $(SRC)/glip_radiometric3.cpp $(INCLUDE)/glip_radiometric3.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_radiometric3.cpp -I$(INCLUDE) -o $(OBJ)/glip_radiometric3.o

$(OBJ)/glip_filtering.o: $(SRC)/glip_filtering.cpp $(INCLUDE)/glip_filtering.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_filtering.cpp -I$(INCLUDE) -o $(OBJ)/glip_filtering.o

$(OBJ)/glip_filtering3.o: $(SRC)/glip_filtering3.cpp $(INCLUDE)/glip_filtering3.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_filtering3.cpp -I$(INCLUDE) -o $(OBJ)/glip_filtering3.o

$(OBJ)/glip_segmentation3.o: $(SRC)/glip_segmentation3.cpp $(INCLUDE)/glip_segmentation3.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_segmentation3.cpp -I$(INCLUDE) -o $(OBJ)/glip_segmentation3.o

$(OBJ)/glip_gradient3.o: $(SRC)/glip_gradient3.cpp $(INCLUDE)/glip_gradient3.h
	$(CXX) $(FLAGS) -c $(SRC)/glip_gradient3.cpp -I$(INCLUDE) -o $(OBJ)/glip_gradient3.o

$(OBJ)/glip_features.o: $(SRC)/glip_features.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_features.cpp -I$(INCLUDE) -o $(OBJ)/glip_features.o

$(OBJ)/glip_attributelist.o: $(SRC)/glip_attributelist.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_attributelist.cpp -I$(INCLUDE) -o $(OBJ)/glip_attributelist.o

$(OBJ)/glip_configfile.o: $(SRC)/glip_configfile.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_configfile.cpp -I$(INCLUDE) -o $(OBJ)/glip_configfile.o

$(OBJ)/glip_geometry.o: $(SRC)/glip_geometry.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_geometry.cpp -I$(INCLUDE) -o $(OBJ)/glip_geometry.o

$(OBJ)/glip_ctree.o: $(SRC)/glip_ctree.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_ctree.cpp -I$(INCLUDE) -o $(OBJ)/glip_ctree.o

$(OBJ)/glip_llw.o: $(SRC)/glip_llw.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_llw.cpp -I$(INCLUDE) -o $(OBJ)/glip_llw.o

$(OBJ)/glip_draw.o: $(SRC)/glip_draw.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_draw.cpp -I$(INCLUDE) -o $(OBJ)/glip_draw.o

$(OBJ)/glip_band.o: $(SRC)/glip_band.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_band.cpp -I$(INCLUDE) -o $(OBJ)/glip_band.o

$(OBJ)/glip_dicom.o: $(SRC)/glip_dicom.cpp
	$(CXX) $(FLAGS) -c $(SRC)/glip_dicom.cpp -I$(INCLUDE) -o $(OBJ)/glip_dicom.o

$(OBJ)/glip_bzlib.o: $(SRC)/glip_bzlib.c
	gcc $(FLAGS) -c $(SRC)/glip_bzlib.c -I$(INCLUDE) -o $(OBJ)/glip_bzlib.o

$(OBJ)/znzlib.o: $(SRC)/znzlib.c
	gcc $(FLAGS) -c $(SRC)/znzlib.c -I$(INCLUDE) -o $(OBJ)/znzlib.o 

$(OBJ)/nifti1_io.o: $(SRC)/nifti1_io.c
	gcc $(FLAGS) -c $(SRC)/nifti1_io.c -I$(INCLUDE) -o $(OBJ)/nifti1_io.o



$(BIN)/areaclosing: $(DEMO)/areaclosing.cpp libglip
	$(CXX) $(FLAGS) $(DEMO)/areaclosing.cpp $(GLIPFLAGS) $(GLIPLIB) -o $(BIN)/areaclosing $(LINKS)

$(BIN)/volumeclosing: $(DEMO)/volumeclosing.cpp libglip
	$(CXX) $(FLAGS) $(DEMO)/volumeclosing.cpp $(GLIPFLAGS) $(GLIPLIB) -o $(BIN)/volumeclosing $(LINKS)

$(BIN)/superpixels: $(DEMO)/superpixels.cpp libglip
	$(CXX) $(FLAGS) $(DEMO)/superpixels.cpp $(GLIPFLAGS) $(GLIPLIB) -o $(BIN)/superpixels $(LINKS)

$(BIN)/supervoxels: $(DEMO)/supervoxels.cpp libglip
	$(CXX) $(FLAGS) $(DEMO)/supervoxels.cpp $(GLIPFLAGS) $(GLIPLIB) -o $(BIN)/supervoxels $(LINKS)

clean:
	rm -f $(LIB)/lib*.a; rm -f $(OBJ)/*.o $(BIN)/*


