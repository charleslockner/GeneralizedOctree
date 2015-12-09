#ifndef __OCTREE_H__
#define __OCTREE_H__

#include <unordered_map>
#include <vector>
#include <Eigen/Dense>

class Cell;

typedef bool(* ObjectCellIntersectionTest)(void * object, Cell * cell);
typedef bool(* ObjectObjectIntersectionTest)(void * objectOut, void * objectIn);

typedef std::vector<void *> ObjectList;
typedef std::vector<Cell *> CellList;
typedef std::unordered_map<void *, CellList> CellMap;
typedef std::pair<void *, CellList> ObjectCellPair;

class Cell {
public:
   Cell(Cell * parent, Eigen::Vector3f lowBound, Eigen::Vector3f highBound);
   ~Cell();

   bool isLeaf();
   void split(ObjectCellIntersectionTest objectInCellTest);

   Eigen::Vector3f lowBound;
   Eigen::Vector3f highBound;
   Eigen::Vector3f center;

   Cell * parent;

   std::vector<Cell *> subcells;
   std::vector<void *> objects;
};

/* Class for efficiently accessing generic objects by location in 3D space.
 */
class Octree {
public:
   Octree(
      Eigen::Vector3f lowBound,
      Eigen::Vector3f highBound,
      unsigned int maxDepth,
      ObjectCellIntersectionTest objectInCellTest
   );
   ~Octree();

   /* Inserts the object to the quadtree and creates nodes as necessary */
   void insert(void * object);

   /* Removes the object from the quadtree */
   void remove(void * object);

   /**
    * Moves the specified object into the correct cells. You must call this if the object is
    * warped or shifted in some way that alters the result of an intersection test.
    */
   void update(void * object);

   /**
    * Removes all data from the octree.
    */
   void clear();

   /**
    *
    */
   void resetWithBounds(Eigen::Vector3f lowBound, Eigen::Vector3f highBound);

   /**
    * Tests for intersection between a specified object already inside the tree, and any other
    * objects within the tree. This is faster than calling testIntersectionOutside.
    * Adds all the objects the specified object collided with to the collisions list parameter.
    * Returns true if there is a collision with another object, false otherwise.
    * If collisions is passed in as NULL, then the octree will not add collision objects to it.
    */
   bool testIntersectionInside(
      void * specObj,
      ObjectObjectIntersectionTest objObjTest,
      ObjectList * collisions
   );

   /**
    * Tests for intersection between a specified object that is not present inside the tree,
    * and any objects within the tree.
    * Adds all the objects the specified object collided with to the collisions list parameter.
    * Returns true if there is a collision with another object, false otherwise.
    * If collisions is passed in as NULL, then the octree will not add collision objects to it.
    */
   bool testIntersectionOutside(
      void * obj,
      ObjectCellIntersectionTest objCellTest,
      ObjectObjectIntersectionTest objObjTest,
      ObjectList * collisions
   );

   /**
    * Test for intersection between a specified object and any other objects within the octree.
    * Adds all the objects the specified object collided with to the collisions list parameter.
    * Returns true if there is a collision with another object, false otherwise.
    * The octree will detect if the specified object is present in the tree. If so, it will
    * already have a reference to all the cells that contain the it, and will traverse only
    * those cells instead of testing all the cells it might be contained in.
    * If objCellTest is passed in as NULL, then Octree will use testIntersectionInside right away.
    * If objCellTest is not NULL and the object is not in the tree, this function calls testIntersectionOutside.
    * If collisions is passed in as NULL, then the octree will not add collision objects to it.
    */
   bool testIntersection(
      void * object,
      ObjectCellIntersectionTest objCellTest,
      ObjectObjectIntersectionTest objObjTest,
      ObjectList * collisions
   );

   Cell * rootCell;

private:
   void addCellToMap(void * object, Cell * cell);
   void insertHelper(void * object, Cell * cell, int lvl);
   void removeSubcellsAndClimbIfEmpty(Cell * cell);
   void clearCellRecursive(Cell * cell);
   bool testIntersectionOutsideHelper(
      void * specObj,
      Cell * cell,
      ObjectCellIntersectionTest objCellTest,
      ObjectObjectIntersectionTest objObjTest,
      ObjectList * collisions
   );

   CellMap cellMap;
   unsigned int maxDepth;
   ObjectCellIntersectionTest objectInCellTest;
};

#endif // __OCTREE_H__