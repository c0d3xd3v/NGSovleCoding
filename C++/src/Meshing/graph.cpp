#include <iostream>
#include "graph.h"

Graph::Graph(int V, GraphCompareCondition *condition)
{
    this->V = V;
    adj.resize(V);
    visited.resize(V, false);
    if(condition != nullptr)
        compareCondition = condition;
    else compareCondition = &defaultCompareCondition;
}

void Graph::addEdge(int v, int w)
{
    // Add w to v’s list.
    adj[v].push_back(w);
}

std::vector<vtkIdType> Graph::BFS(int s)
{
    std::vector<vtkIdType> ret;
    // Mark all the vertices as not visited

    // Create a queue for BFS
    std::list<int> queue;

    // Mark the current node as visited and enqueue it
    visited[s] = true;
    queue.push_back(s);

    while (!queue.empty()) {

        // Dequeue a vertex from queue and print it
        int k = queue.front();
        //std::cout << k << " ";
        ret.push_back(k);
        queue.pop_front();

        // Get all adjacent vertices of the dequeued
        // vertex s.
        // If an adjacent has not been visited,
        // then mark it visited and enqueue it
        for (auto adjacent : adj[k]) {
            if (!visited[adjacent]) {
                visited[adjacent] = true;
                if(compareCondition->compare(k, adjacent))
                    queue.push_back(adjacent);
            }
        }
    }

    return ret;
}

void Graph::reset()
{
    //visited.resize(V, false);
    for(int i = 0; i < visited.size(); i++)
        visited[i] = false;
}
