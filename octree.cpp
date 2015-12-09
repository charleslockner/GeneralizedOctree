
#include "octree.h"

#include "matrix_math.h"

#include <algorithm>

/* subcell order:
 * 0: (-,-,-)
 * 1: (-,-,+)
 * 2: (-,+,-)
 * 3: (-,+,+)
 * 4: (+,-,-)
 * 5: (+,-,+)
 * 6: (+,+,-)
 * 7: (+,+,+)
 */

Cell::Cell(Cell * parent, Eigen::Vector3f lowBound, Eigen::Vector3f highBound) {
   this->parent = parent;
   this->lowBound = lowBound;
   this->highBound = highBound;
   this->center = Eigen::Vector3f(
      (lowBound(0) + highBound(0)) / 2.0f,
      (lowBound(1) + highBound(1)) / 2.0f,
      (lowBound(2) + highBound(2)) / 2.0f
   );
   this->subcells = std::vector<Cell *>();
   this->objects = std::vector<void *>();
}

Cell::~Cell() {}

bool Cell::isLeaf() {
   return subcells.size() == 0;
}

// Subdivide cell and move all objects into created subcells
void Cell::split(ObjectCellIntersectionTest objectInCellTest) {
   /* Create each octant/subcell */
   // Octant: (-,-,-) <=> (0,0,0)
   subcells.push_back(new Cell(
      this,
      lowBound,
      center
   ));
   // Octant: (-,-,0) <=> (0,0,+)
   subcells.push_back(new Cell(
      this,
      Eigen::Vector3f(lowBound(0), lowBound(1), center(2)),
      Eigen::Vector3f(center(0), center(1), highBound(2))
   ));
   // Octant: (-,0,-) <=> (0,+,0)
   subcells.push_back(new Cell(
      this,
      Eigen::Vector3f(lowBound(0), center(1), lowBound(2)),
      Eigen::Vector3f(center(0), highBound(1), center(2))
   ));
   // Octant: (-,0,0) <=> (0,+,+)
   subcells.push_back(new Cell(
      this,
      Eigen::Vector3f(lowBound(0), center(1), center(2)),
      Eigen::Vector3f(center(0), highBound(1), highBound(2))
   ));
   // Octant: (0,-,-) <=> (+,0,0)
   subcells.push_back(new Cell(
      this,
      Eigen::Vector3f(center(0), lowBound(1), lowBound(2)),
      Eigen::Vector3f(highBound(0), center(1), center(2))
   ));
   // Octant: (0,-,0) <=> (+,0,+)
   subcells.push_back(new Cell(
      this,
      Eigen::Vector3f(center(0), lowBound(1), center(2)),
      Eigen::Vector3f(highBound(0), center(1), highBound(2))
   ));
   // Octant: (0,0,-) <=> (+,+,0)
   subcells.push_back(new Cell(
      this,
      Eigen::Vector3f(center(0), center(1), lowBound(2)),
      Eigen::Vector3f(highBound(0), highBound(1), center(2))
   ));
   // Octant: (0,0,0) <=> (+,+,+)
   subcells.push_back(new Cell(
      this,
      center,
      highBound
   ));
}

Octree::Octree(
   Eigen::Vector3f lowBound,
   Eigen::Vector3f highBound,
   unsigned int maxDepth,
   ObjectCellIntersectionTest objectInCellTest
) {
   rootCell = new Cell(NULL, lowBound, highBound);
   this->maxDepth = maxDepth;
   this->objectInCellTest = objectInCellTest;
   cellMap = CellMap();
}

Octree::~Octree() {

}


// Add the specified cell to the object's list of cell references in cellMap
void Octree::addCellToMap(void * object, Cell * cell) {
   CellMap::iterator it = cellMap.find(object);
   if (it == cellMap.end()) {
      // printf("adding cell to map for first time, cell = %p\n", cell);
      CellList cells = CellList();
      cells.push_back(cell);
      cellMap.insert(ObjectCellPair(object, cells));
   } else {
      // printf("adding cell to map (not first time), cell = %p\n", cell);
      it->second.push_back(cell);
   }
}

// Recursive helper function that adds the object to each leafcell that will contain it.
// Also splits any leaf cell that will contain the object and will have more than max objects in it
void Octree::insertHelper(void * object, Cell * cell, int lvl) {
   if (objectInCellTest(object, cell)) {
      if (cell->isLeaf()) {   // Is a leaf cell
         if (lvl == maxDepth) { // This cell is at max depth
            // So, add the object to the cell and be done with this recursion
            cell->objects.push_back(object);
            addCellToMap(object, cell);
         } else { // this cell can be split since not at max depth
            // So, split the cell and recurse one more time
            cell->split(objectInCellTest);
            for (int i = 0; i < 8; i++) {
               insertHelper(object, cell->subcells[i], lvl+1);
            }
         }
      } else {                // Is not a leaf cell
         // So, keep recursing
         for (int i = 0; i < 8; i++) {
            insertHelper(object, cell->subcells[i], lvl+1);
         }
      }
   }
}

void Octree::insert(void * object) {
   insertHelper(object, rootCell, 0);
}

void Octree::removeSubcellsAndClimbIfEmpty(Cell * cell) {
   // Check to see if we can delete further cells higher up
   for (int i = 0; i < 8; i++) {
      if (!cell->subcells[i]->isLeaf() || cell->subcells[i]->objects.size() > 0) {
         return ;
      }
   }

   // Will only get here if no subcells contain an object and all subcells are not leaves
   for (int i = 0; i < 8; i++) {
      delete(cell->subcells[i]);
   }
   cell->subcells.clear();

   // Keep recursing if we're not at the root yet
   if (cell->parent != NULL) {
      removeSubcellsAndClimbIfEmpty(cell->parent);
   }
}

void Octree::remove(void * specObj) {
   CellMap::iterator cellIt = cellMap.find(specObj);
   if (cellIt == cellMap.end()) {
      fprintf(stderr, "Octree::removeObjectFromCells WARNING: the specified object was not found in the tree.\n");
      return ;
   }

   CellList& cells = cellIt->second;
   int numCells = cells.size();
   for (int i = 0; i < numCells; i++) {
      Cell * cell = cells[i];

      ObjectList& objs = cell->objects;
      objs.erase(std::remove(objs.begin(), objs.end(), specObj), objs.end());

      if (cell->parent) {
         removeSubcellsAndClimbIfEmpty(cell->parent);
      }
   }

   cellMap.erase(specObj);
}

void Octree::update(void * specObj) {
   remove(specObj);
   insert(specObj);
}

void Octree::clearCellRecursive(Cell * cell) {
   if (!cell->isLeaf()) {
      for (int i = 0; i < 8; i++) {
         clearCellRecursive(cell->subcells[i]);
      }
   }
   if (cell != rootCell) {
      delete(cell);
   }
}

void Octree::clear() {
   clearCellRecursive(rootCell);
   rootCell->subcells.clear();
   rootCell->objects.clear();
   cellMap.clear();
}

void Octree::resetWithBounds(Eigen::Vector3f lowBound, Eigen::Vector3f highBound) {
   clearCellRecursive(rootCell);
   rootCell->subcells.clear();
   rootCell->objects.clear();

   rootCell->lowBound = lowBound;
   rootCell->highBound = highBound;
   rootCell->center = Eigen::Vector3f(
      (lowBound(0) + highBound(0)) / 2.0f,
      (lowBound(1) + highBound(1)) / 2.0f,
      (lowBound(2) + highBound(2)) / 2.0f
   );

   for (CellMap::iterator it = cellMap.begin(); it != cellMap.end(); it++) {
      it->second.clear();
      insert(it->first);
   }
}

bool Octree::testIntersectionOutsideHelper(
   void * specObj,
   Cell * cell,
   ObjectCellIntersectionTest objCellTest,
   ObjectObjectIntersectionTest objObjTest,
   ObjectList * collisions
) {
   bool hasCollision = false;

   if (objCellTest(specObj, rootCell)) {
      if (cell->isLeaf()) {
         int numObjects = cell->objects.size();
         for (int i = 0; i < numObjects; i++) {
            void * obj = cell->objects[i];
            if (obj != specObj && objObjTest(specObj, obj)) {
               if (collisions != NULL) {
                  collisions->push_back(obj);
               }
               hasCollision = true;
            }
         }
      } else {
         for (int i = 0; i < 8; i++) {
            hasCollision |= testIntersectionOutsideHelper(specObj, cell->subcells[i], objCellTest, objObjTest, collisions);
         }
      }
   }

   return hasCollision;
}

bool Octree::testIntersectionInside(
   void * specObj,
   ObjectObjectIntersectionTest objObjTest,
   ObjectList * collisions
) {
   CellMap::iterator it = cellMap.find(specObj);
   if (it == cellMap.end()) {
      fprintf(stderr, "Octree::testIntersectionInside WARNING: the specified object was not found in the tree.\n");
      return false;
   }

   CellList& cells = it->second;
   bool hasCollision = false;

   int numCells = cells.size();
   for (int i = 0; i < numCells; i++) {
      Cell * cell = cells[i];
      int numObjs = cell->objects.size();
      for (int j = 0; j < numObjs; j++) {
         void * obj = cell->objects[j];
         if (obj != specObj && objObjTest(specObj, obj)) {
            if (collisions != NULL) {
               collisions->push_back(obj);
            }
            hasCollision = true;
         }
      }
   }

   return hasCollision;
}

bool Octree::testIntersectionOutside(
   void * obj,
   ObjectCellIntersectionTest objCellTest,
   ObjectObjectIntersectionTest objObjTest,
   ObjectList * collisions
) {
   return testIntersectionOutsideHelper(obj, rootCell, objCellTest, objObjTest, collisions);
}

/**
 * Test for intersection between a specified object and any other objects within the octree.
 * Adds all the objects the specified object collided with to the collisions list parameter.
 * Returns true if there is a collision with another object, false otherwise.
 * The octree will detect if the specified object is present in the tree. If so, it will
 * already have a reference to all the cells that contain the it, and will traverse only
 * those cells instead of testing all the cells it might be contained in.
 * If objCellTest is passed in as NULL, then Octree will use testIntersectionInside right away.
 */
bool Octree::testIntersection(
   void * object,
   ObjectCellIntersectionTest objCellTest,
   ObjectObjectIntersectionTest objObjTest,
   ObjectList * collisions
) {
   if (objCellTest == NULL || cellMap.find(object) != cellMap.end()) {
      return testIntersectionInside(object, objObjTest, collisions);
   } else {
      return testIntersectionOutside(object, objCellTest, objObjTest, collisions);
   }
}
