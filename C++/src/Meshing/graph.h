#ifndef GRAPH_H
#define GRAPH_H

#include <list>
#include <vector>
#include <vtkType.h>

class GraphCompareCondition {
public:
    virtual bool compare(vtkIdType f1, vtkIdType f2) = 0;
};

class DefaultGraphCompareCondition : public GraphCompareCondition {
public:
    bool compare(vtkIdType f1, vtkIdType f2) {
        return true;
    }
};

// This class represents a directed graph using
// adjacency list representation
class Graph {
private:
    // No. of vertices
    int V;
    std::vector<bool> visited;
    // Pointer to an array containing adjacency lists
    std::vector<std::list<vtkIdType> > adj;
    DefaultGraphCompareCondition defaultCompareCondition;
    GraphCompareCondition *compareCondition;

public:
    // Constructor
    Graph(int V, GraphCompareCondition *condition  = nullptr);

    // Function to add an edge to graph
    void addEdge(int v, int w);

    // Prints BFS traversal from a given source s
    std::vector<vtkIdType> BFS(int s);

    void reset();
};

#endif // GRAPH_H
