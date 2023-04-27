import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

public class Graph {

    class Node {

        //Vertex Class Fields
        private List<Edge> edgeList;
        private final String id;
        private boolean encountered;

        //Vertex Class Constructor
        public Node(String id){
            this.id = id;
            this.edgeList = new ArrayList<>();
            encountered = false;
        }
    }

    class Edge {

        //Edge Class Fields
        private List<Node> edgeNodes;
        private Node startNode;
        private Node endNode;

        //Edge Class Constructor
        public Edge(Node startNode, Node endNode) {
            this.edgeNodes = new ArrayList<>();
            this.startNode = startNode;
            this.endNode = endNode;
            edgeNodes.add(startNode);
            edgeNodes.add(endNode);
        }
    }

    //Graph Class Fields
    private List<Node> adjList;
    private int totalNumEdges;

    //Graph Class Constructor
    public Graph(){
            this.adjList = new ArrayList<>();
            totalNumEdges = 0;
        }

        //Method to add node to graph
        public boolean addNode(String nodeName){
            for (Node node : adjList)
            if (node.id == nodeName) {
                System.out.println("\nNode " + nodeName + " is already in the graph!");
                return false;
            }

            Node addedNode = new Node(nodeName);
            adjList.add(addedNode);
            return true;
        }

        //Helper method to remove duplicates in a String
        public String[] removeDuplicates(String[]  input){
            input = new HashSet<String>(Arrays.asList(input)).toArray(new String[0]); //Used built in java hash set because it removes duplicates automatically
            return input;
  }

        //Method to add multiple nodes into graph where every String in the input array will added as a node
        public boolean addNodes(String[] nodeNames) {
            boolean alrInGraph = false;
            String[] noDupNodeNames = removeDuplicates(nodeNames);
            List<String> noDupNodeNames2 = new ArrayList<>();
            for (String s : noDupNodeNames)
                noDupNodeNames2.add(s);

            for (String nodeName : noDupNodeNames) {
                for (Node node2 : adjList) {
                    if (node2.id == nodeName) {
                        System.out.println("\nNode " + nodeName + " is already in the graph!");
                        alrInGraph = true;
                        noDupNodeNames2.remove(node2.id);
                    }
                }

                }
            for(String s : noDupNodeNames2) {
                Node addedNode = new Node(s);
                adjList.add(addedNode);
            }

            //If any of the input nodes are already in the graph, we return false
            if (alrInGraph == true)
                return false;

            //We return true only if all input nodes were not already in the graph and were successfully added to the graph
            return true;
        }

        //Helper Method to fetch Node with the same name as the string input
        public Node fetchNode(String input){
            int nodeIndex = 0;


            //After this for loop is run, adjList.get(startNodeIndex) is the startNode
            for(int i = 0; i < adjList.size(); i++){
                if (input == adjList.get(i).id){
                    nodeIndex = i;
                    break;
                }
            }

            Node foundNode = adjList.get(nodeIndex);
            return foundNode;
        }

        //Method to add and edge that connects the startNode and endNode
        public boolean addEdge(String startNode, String endNode) {

            //Checking if startNode exists in the graph
            boolean startNodeIncluded = false;
            boolean endNodeIncluded = false;
            for (Node node : adjList)
                if (node.id == startNode)
                    startNodeIncluded = true;
                if (startNodeIncluded == false) {
                    System.out.println("\nCannot create an edge including node " + startNode + " because " + startNode + " doesn't exist");
                    return false;
                }

                //Checking if endNode exists in the graph
            for (Node node : adjList)
                if (node.id == endNode)
                    endNodeIncluded = true;
            if (endNodeIncluded == false) {
                System.out.println("\nCannot create an edge including node " + endNode + " because " + endNode + " doesn't exist");
                return false;
            }

            //Now we know node exists, fetch the node with the name of the string input
            Node newStartNode = fetchNode(startNode);
            Node newEndNode = fetchNode(endNode);

            //Create new edge containing the two nodes
            Edge newEdge = new Edge(newStartNode, newEndNode);

            //Add the edge to the edge list of the start node given in the input
            newStartNode.edgeList.add(newEdge);
            System.out.println("\nEdge ( " + startNode + "<--->" + endNode + " ) has been added!");
            totalNumEdges++;


            //This if statement makes sure that if startNode and endNode are the same (self-loop), the edge won't be added to the same vertex twice
            if (startNode != endNode) {

                //Creating the edge that will be placed in the edgeList of the endNode
                Edge newEdge2 = new Edge(newEndNode, newStartNode);

                //Add the edge to the edge list of the start node given in the input
                newEndNode.edgeList.add(newEdge2);

            }
            return true;
        }

        //Method that adds multiple edges from startNode to every endNode in a String array
        public boolean addEdges(String startNode, String[] endNodes){

        //Checking if startNode exists in the graph
            boolean startNodeIncluded = false;
            for (Node node : adjList)
                if (node.id == startNode)
                    startNodeIncluded = true;
            if (startNodeIncluded == false) {
                System.out.println("\nCannot create an edge including node " + startNode + " because " + startNode + " doesn't exist");
                return false;
            }

            //Checking if any endNodes exist in the graph
            List<String> includedEndNodes = new ArrayList<>(); //Creating arraylist to keep track of nodes that exist
            for (String endNode : endNodes)
                includedEndNodes.add(endNode);

            int numEndNodeIncluded = 0;
            for (String endNode : endNodes) {
                boolean isEndNodeIncluded = false;
                for (Node node : adjList){
                    if (node.id == endNode) {
                        isEndNodeIncluded = true;
                        numEndNodeIncluded++;
                    }
                    }
                if (isEndNodeIncluded == false) {
                    includedEndNodes.remove(endNode);
                    System.out.println("\nCannot create an edge including node " + endNode + " because " + endNode + " doesn't exist");

                }
            }
            if (numEndNodeIncluded == 0){
                System.out.println("\nAll end nodes don't exist, so no edges were created");
                return false;
            }

            //All edges will be created in this for loop
            for (String endNode : includedEndNodes) {

                //Create new nodes corresponding to inputs
                Node newStartNode = fetchNode(startNode);
                Node newEndNode = fetchNode(endNode);

                //Create new edge containing the two nodes
                Edge newEdge = new Edge(newStartNode, newEndNode);

                //Add the edge to the edge list of the start node given in the input
                newStartNode.edgeList.add(newEdge);
                System.out.println("\nEdge ( " + startNode + "<--->" + endNode + " ) has been added!");
                totalNumEdges++;


                //This if statement makes sure that if startNode and endNode are the same (self-loop), the edge won't be added to the same vertex twice
                if (endNode != startNode) {

                    //Creating the edge that will be placed in the edgeList of the endNode
                    Edge newEdge2 = new Edge(newEndNode, newStartNode);

                    //Add the edge to the edge list of the start node given in the input
                    newEndNode.edgeList.add(newEdge2);

                }
            }
            return true;
        }

        //Helper Method to delete the other end of an edge in a nodes edgeList
        public void deleteOtherEdgeEnd(String node) {
            Node initialNode = fetchNode(node);
            for (Edge edge : initialNode.edgeList) {
                Node nodeWithOtherEdge = fetchNode(edge.endNode.id);
                for (Edge edge2 : nodeWithOtherEdge.edgeList) {
                    if (edge2.endNode.id == node) {
                        nodeWithOtherEdge.edgeList.remove(edge2);
                        totalNumEdges--;
                        break;
                    }
                }
            }
        }

        //Method to remove a node and all of its edges from the graph
        public boolean removeNode(String nodeName) {

            for (Node n : adjList) {
                if (n.id == nodeName) {
                    Node node = fetchNode(nodeName);
                    deleteOtherEdgeEnd(nodeName);
                    adjList.remove(node);
                    System.out.println("\nNode " + node.id + " has been removed!");
                    return true;
                }
            }

            System.out.println("\nNode " + nodeName + " cannot be removed because node " + nodeName + " doesn't exist");
            return false;
        }

        //Method to remove multiple nodes from the graph
        public boolean removeNodes(String[] nodeNames) {
            String[] noDupNodeNames = removeDuplicates(nodeNames);
            List<String> noDupNodeNames2 = new ArrayList<>();
            for (String s : noDupNodeNames)
                noDupNodeNames2.add(s);
            int initialInputSize = noDupNodeNames2.size();
            int initialAdjListSize = adjList.size();

            //Checking if any of the nodeNames are already in the graph
            for (int i = 0; i < initialInputSize; i++) {
                String nodeName = noDupNodeNames2.get(i);
                boolean alrInGraph = false;
                for (Node node : adjList) {
                    if (node.id == nodeName) {
                        alrInGraph = true;
                        break;
                    }
                }
                if (alrInGraph == false) {
                    noDupNodeNames2.remove(i);
                    i--;
                    System.out.println("\nNode " + nodeName + " cannot be removed because node " + nodeName + " doesn't exist");
                }
                initialInputSize = noDupNodeNames2.size();
            }

            //If all input nodes are not in the graph, we return false
            if (noDupNodeNames2.size() == 0) {
                System.out.println("No nodes were removed");
                return false;
            }

            //the allNodesRemoved array will be used to print all the nodes we removed when method is finished
            String[] allNodesRemoved = new String[noDupNodeNames2.size()];
            for (int i = 0; i < noDupNodeNames2.size(); i++) {
                allNodesRemoved[i] = noDupNodeNames2.get(i);
            }

            //if adj list does contain nodeName...
            List<String> nodesRemoved = new ArrayList<>();
            for (String nodeName : noDupNodeNames2) {
                Node node = fetchNode(nodeName);
                deleteOtherEdgeEnd(nodeName);
                adjList.remove(node);
                nodesRemoved.add(nodeName);
            }

            //If all nodes were removed
            if (adjList.size() == initialAdjListSize - noDupNodeNames2.size()) {
                System.out.println("\nAll nodes " + Arrays.toString(allNodesRemoved) + " were removed!");
                return true;
            }

            //If some, but not all, nodes were removed
            String[] printNodesRemoved = new String[nodesRemoved.size()];
            for (int i = 0; i < nodesRemoved.size(); i++)
                printNodesRemoved[i] = nodesRemoved.get(i);
            System.out.println("\nOnly nodes " + printNodesRemoved + " were removed!");
            return false;
        }

            //Method to print graph in adjacency list format
            public void printGraph() {
                System.out.println("\n");
                System.out.println("                         Graph Adjacency List");
                System.out.println("----------------------------------------------------------------------");
                for (int i = 0; i < adjList.size(); i++) {
                    System.out.print("Vertex: " + adjList.get(i).id + "     Edges: ");
                    for (int j = 0; j < adjList.get(i).edgeList.size(); j++) {
                        System.out.print(adjList.get(i).edgeList.get(j).endNode.id + " || ");
                    }
                    System.out.println("\n----------------------------------------------------------------------");
                }
                System.out.println("----------------------------------------------------------------------");
                System.out.println("Total Number of Vertices: " + adjList.size());
                System.out.println("Total Number of Edges: " + totalNumEdges);
            }

            //----------------------------------------------------------------------------------------

        //Search Algorithms!


        /*
        Helper method to find neighbor node that comes earliest in alphabetical order
            -Only works if input node has neighbors
         */

        public String[]  orderNeighborsAlphaOrder(String input) {
            Node firstNode = fetchNode(input);
            List<String> neighborNodeNames = new ArrayList<>();
            boolean unvisitedNodeExists = false;
            String[] inputArr = new String[]{"input"};

            //For loop to find the name of all the input node's neighbor nodes that haven't bee visited

                for (Edge edge : firstNode.edgeList) {
                    if (edge.endNode.encountered == false) {
                        unvisitedNodeExists = true;
                        String s = edge.endNode.id.replaceAll("[^a-zA-Z]", "").toLowerCase();
                        neighborNodeNames.add(s);
                    }
                }
                if (unvisitedNodeExists == false)
                    return inputArr;

                String[] neighborNodeNamesArr = arrListToArray(neighborNodeNames);
                neighborNodeNamesArr = Stream.of(neighborNodeNamesArr).sorted(Comparator.reverseOrder()).toArray(String[]::new);
                return neighborNodeNamesArr;


        }

        /*
        Helper method to find neighbor node that comes earliest in reverse alphabetical order
            -Only works if input node has neighbors
         */

    public String[]  orderNeighborsReverseAlphaOrder(String input) {
        Node firstNode = fetchNode(input);
        List<String> neighborNodeNames = new ArrayList<>();
        boolean unvisitedNodeExists = false;
        String[] inputArr = new String[]{"input"};

        //For loop to find the name of all the input node's neighbor nodes that haven't bee visited

        for (Edge edge : firstNode.edgeList) {
            if (edge.endNode.encountered == false) {
                unvisitedNodeExists = true;
                String s = edge.endNode.id.replaceAll("[^a-zA-Z]", "").toLowerCase();
                neighborNodeNames.add(s);
            }
        }
        if (unvisitedNodeExists == false)
            return inputArr;

        String[] neighborNodeNamesArr = arrListToArray(neighborNodeNames);
        neighborNodeNamesArr = Stream.of(neighborNodeNamesArr).sorted().toArray(String[]::new);
        return neighborNodeNamesArr;
    }

        //Helper method to copy array list to an array
        public String[] arrListToArray(List<String> s){
            String[] arr = new String[s.size()];
            for (int i = 0; i < s.size(); i++)
                arr[i] = s.get(i);
            return arr;
        }


        /*
        Method for Depth First Search Algorithm

        -This method returns the path found using the algorithm between node "from" and node "to"

        -The neighborOrder can be either "alphabetical" or "reverse".
            -alphabetical means that the order in which neighbors are considered is by alphabetical order
            -reverse means that the order in which neighbors are considered is by REVERSE alphabetical order
         */
        public String[] DFS(String from, String to, String neighborOrder) {
            Stack<Node> stack = new Stack<Node>();
            String[] emptyArr = new String[0];

            //Checking if node "from" exists
            boolean nodeFromExists = false;
            for (Node node : adjList){
                if (node.id == from)
                    nodeFromExists = true;
            }
            if (nodeFromExists == false) {
                System.out.println("Start node " + from + " doesn't exist");
                return emptyArr;
            }

            //Checking if node "to" exists
            boolean nodeToExists = false;
            for (Node node : adjList){
                if (node.id == to)
                    nodeToExists = true;
            }
            if (nodeToExists == false) {
                System.out.println("Target node " + to + " doesn't exist");
                return emptyArr;
            }

            //Checking if neighbor order is "alphabetical" or "reverse"
            if (neighborOrder != "alphabetical" & neighborOrder != "reverse"){
                System.out.println("Neighbor node search order needs to be alphabetical or reverse");
                return emptyArr;
            }

            //Now we know that all inputs are valid, we can begin DFS
            List<String> visited = new ArrayList<>();
            Node startNode = fetchNode(from);
            stack.push(startNode);

            if (neighborOrder == "alphabetical") {
                while (!stack.isEmpty()) {
                    Node current = stack.pop();
                    if (current.id == to) {
                        System.out.println("\nThe path from node " + from + " to node " + to + " is " + Arrays.toString(arrListToArray(visited)) + " via the Depth First Search algorithm");
                        return arrListToArray(visited);
                    } else {
                        if (!visited.contains(current.id)) {
                            visited.add(current.id);
                            /* This String array we're iterating through is actually reverse alphabetical order because we want the first string alphabetically
                            to be pushed last. The elements at the top of the stack are the ones that we check first
                             */
                            for (String neighbor : orderNeighborsAlphaOrder(current.id)) {
                                stack.push(fetchNode(neighbor));
                            }
                        }
                    }
                }
            }

            //if neighborOrder input isn't alphabetical, it has to be reverse order
            while (!stack.isEmpty()) {
                Node current = stack.pop();
                if (current.id == to) {
                    visited.add(current.id);
                    System.out.println("\nThe path from node " + from + " to node " + to + " is " + Arrays.toString(arrListToArray(visited)) + " via the Depth First Search algorithm");
                    return arrListToArray(visited);
                } else {
                    if (!visited.contains(current.id)) {
                        visited.add(current.id);
                        for (String neighbor : orderNeighborsReverseAlphaOrder(current.id)) {
                            stack.push(fetchNode(neighbor));
                        }
                    }
                }
            }

            System.out.println("\nPath doesn't exist");
            return emptyArr;
        }

        public String[] BFS(String from, String to, String neighborOrder) {
            Queue<Node> queue = new LinkedList<>();
            List<String> visited = new ArrayList<>();
            String[] emptyArr = new String[0];


            //Checking if node "from" exists
            boolean nodeFromExists = false;
            for (Node node : adjList) {
                if (node.id == from)
                    nodeFromExists = true;
            }
            if (nodeFromExists == false) {
                System.out.println("Start node " + from + " doesn't exist");
                return emptyArr;
            }

            //Checking if node "to" exists
            boolean nodeToExists = false;
            for (Node node : adjList) {
                if (node.id == to)
                    nodeToExists = true;
            }
            if (nodeToExists == false) {
                System.out.println("Target node " + to + " doesn't exist");
                return emptyArr;
            }

            //Checking if neighbor order is "alphabetical" or "reverse"
            if (neighborOrder != "alphabetical" & neighborOrder != "reverse") {
                System.out.println("Neighbor node search order needs to be alphabetical or reverse");
                return emptyArr;
            }

            Node current = fetchNode(from);
            queue.add(current);
            visited.add(from);

            while (!queue.isEmpty()) {
                current = queue.poll();

                if (neighborOrder == "alphabetical") {
                    for (String neighbor : orderNeighborsReverseAlphaOrder(current.id)) {
                        if (!visited.contains(neighbor)) {
                            queue.add(fetchNode(neighbor));
                            visited.add(neighbor);
                        }

                        if (neighbor == to) {
                            System.out.println("\nThe path from node " + from + " to node " + to + " is " + Arrays.toString(arrListToArray(visited)) + " via the Depth First Search algorithm");
                            return arrListToArray(visited);
                        }
                    }
                }

                for (String neighbor : orderNeighborsAlphaOrder(current.id)) {
                    if (!visited.contains(neighbor)) {
                        queue.add(fetchNode(neighbor));
                        visited.add(neighbor);
                    }

                    if (neighbor == to) {
                        System.out.println("\nThe path from node " + from + " to node " + to + " is " + Arrays.toString(arrListToArray(visited)) + " via the Depth First Search algorithm");
                        return arrListToArray(visited);
                    }
                }


            }

            System.out.println("\nPath doesn't exist");
            return emptyArr;
        }

        //Method to return shortest path between node from and node to
        public String[] shortestPath(String from, String to) {
            String[] emptyArr = new String[0];

            //Checking if node "from" exists
            boolean nodeFromExists = false;
            for (Node node : adjList) {
                if (node.id == from)
                    nodeFromExists = true;
            }
            if (nodeFromExists == false) {
                System.out.println("Start node " + from + " doesn't exist");
                return emptyArr;
            }

            //Checking if node "to" exists
            boolean nodeToExists = false;
            for (Node node : adjList) {
                if (node.id == to)
                    nodeToExists = true;
            }
            if (nodeToExists == false) {
                System.out.println("Target node " + to + " doesn't exist");
                return emptyArr;
            }

            // Initialize distances and predecessors
            Map<String, Integer> distances = new HashMap<>();
            Map<String, String> predecessors = new HashMap<>();
            for (Node node : adjList) {
                distances.put(node.id, Integer.MAX_VALUE);
                predecessors.put(node.id, null);
            }
            distances.put(from, 0);

            // Initialize priority queue for Dijkstra's algorithm
            PriorityQueue<Map.Entry<String, Integer>> queue = new PriorityQueue<>((x, y) -> distances.get(x.getKey()) - distances.get(y.getKey()));
            queue.add(new AbstractMap.SimpleEntry<>(from, 0));

            // Dijkstra's algorithm
            while (!queue.isEmpty()) {
                Map.Entry<String, Integer> current = queue.poll();
                String node = current.getKey();
                int distance = current.getValue();

                if (distances.get(node) < distance)
                    continue;

                for (Edge edge : fetchNode(current.getKey()).edgeList) {
                    String neighbor = edge.endNode.id;
                    int newDistance = distances.get(node) + 1;

                    if (newDistance < distances.get(neighbor)) {
                        distances.put(neighbor, newDistance);
                        predecessors.put(neighbor, node);
                        queue.add(new AbstractMap.SimpleEntry<>(neighbor, newDistance));
                    }
                }
            }

            // Backtrack to construct the shortest path
            List<String> path = new ArrayList<>();
            String current = to;
            while (current != null) {
                path.add(current);
                current = predecessors.get(current);
            }
            Collections.reverse(path);
            String[] output = path.toArray(new String[path.size()]);

            System.out.println("\nThe shortest path from node " + from + " to node " + to + " is " + Arrays.toString(output) + " via Dijkstra's algorithm");
            return output;
        }

    //fetchNode helper method, except it takes in a specific adjacency list as in input
    public Node fetchNode2 (String input, List<Node> adjList){
        int nodeIndex = 0;


        //After this for loop is run, adjList.get(startNodeIndex) is the startNode
        for(int i = 0; i < adjList.size(); i++){
            if (input == adjList.get(i).id){
                nodeIndex = i;
                break;
            }
        }

        Node foundNode = adjList.get(nodeIndex);
        return foundNode;
    }

    //Method to return second shortest path between node from and node to
    public String[] secondShortestPath(String from, String to){
        String[] emptyArr = new String[0];

        //Checking if node "from" exists
        boolean nodeFromExists = false;
        for (Node node : adjList) {
            if (node.id == from)
                nodeFromExists = true;
        }
        if (nodeFromExists == false) {
            System.out.println("Start node " + from + " doesn't exist");
            return emptyArr;
        }

        //Checking if node "to" exists
        boolean nodeToExists = false;
        for (Node node : adjList) {
            if (node.id == to)
                nodeToExists = true;
        }
        if (nodeToExists == false) {
            System.out.println("Target node " + to + " doesn't exist");
            return emptyArr;
        }

            String[] firstShortestPath = shortestPath(from, to);
            List<Node> adjList2 = adjList;

        //This for loop deletes the edge between input node "to" and the node that came before it in the path returned by the shortest algorithm method
        for (int i = 0; i < adjList2.size(); i++){
            if (adjList2.get(i).id == firstShortestPath[firstShortestPath.length-1]) {
                for (int j = 0; j < adjList2.get(i).edgeList.size(); j++) {
                    if (adjList2.get(i).edgeList.get(j).endNode.id == firstShortestPath[firstShortestPath.length - 2])
                        adjList2.get(i).edgeList.remove(j);
                }
            }

            if (adjList2.get(i).id == firstShortestPath[firstShortestPath.length-2]) {
                for (int j = 0; j < adjList2.get(i).edgeList.size(); j++) {
                    if (adjList2.get(i).edgeList.get(j).endNode.id == firstShortestPath[firstShortestPath.length - 1])
                        adjList2.get(i).edgeList.remove(j);
                }
            }
        }


        /*Now we run the same the path finding algorithm we used in the shortestPath method, except instead of searching through the original adjList, it searches through
        the adjList2 which is just adjList without the final edge in the first shortest path
         */

        // Initialize distances and predecessors
        Map<String, Integer> distances = new HashMap<>();
        Map<String, String> predecessors = new HashMap<>();
        for (Node node : adjList2) {
            distances.put(node.id, Integer.MAX_VALUE);
            predecessors.put(node.id, null);
        }
        distances.put(from, 0);

        // Initialize priority queue for Dijkstra's algorithm
        PriorityQueue<Map.Entry<String, Integer>> queue = new PriorityQueue<>((x, y) -> distances.get(x.getKey()) - distances.get(y.getKey()));
        queue.add(new AbstractMap.SimpleEntry<>(from, 0));

        // Dijkstra's algorithm
        while (!queue.isEmpty()) {
            Map.Entry<String, Integer> current = queue.poll();
            String node = current.getKey();
            int distance = current.getValue();

            if (distances.get(node) < distance)
                continue;

            for (Edge edge : fetchNode2(current.getKey(), adjList2).edgeList) {
                String neighbor = edge.endNode.id;
                int newDistance = distances.get(node) + 1;

                if (newDistance < distances.get(neighbor)) {
                    distances.put(neighbor, newDistance);
                    predecessors.put(neighbor, node);
                    queue.add(new AbstractMap.SimpleEntry<>(neighbor, newDistance));
                }
            }
        }

        // Backtrack to construct the shortest path
        List<String> path = new ArrayList<>();
        String current = to;
        while (current != null) {
            path.add(current);
            current = predecessors.get(current);
        }
        Collections.reverse(path);
        String[] output = path.toArray(new String[path.size()]);

        System.out.println("\nThe second shortest path from node " + from + " to node " + to + " is " + Arrays.toString(output) + " via Dijkstra's algorithm");
        return output;
    }


        public static void main(String[] args){

        Graph mufu = new Graph();
        String[] test = new String[]{"mufu", "ayo", "breana", "ope", "josh", "kinny", "will", "ally", "faith"};
        String[] test2 = new String[]{"ayo", "breana", "ope"};
        String[] test3 = new String[]{"josh", "kinny"};
        String[] test4 = new String[]{"will", "all of case"};
        mufu.addNodes(test);
        mufu.addEdges("mufu", test2);
        mufu.addEdges("ayo", test3);
        mufu.addEdges("breana", test4);
        mufu.addEdge("ope", "faith");
        mufu.addNode("maxine");
        mufu.addEdge("maxine", "kinny");
        mufu.addEdge("maxine", "mufu");
        mufu.printGraph();
        mufu.shortestPath("mufu", "maxine");
        mufu.secondShortestPath("mufu", "maxine");
        mufu.BFS("mufu", "will", "reverse");

        }
}







