#ifndef QGSTACK_H
#define QGSTACK_H

#include "qgstring.hh"
#include "assert.h"

class qgstack
{
  public:
   qgstack(); //constructor
   ~qgstack(); // destructor
   void put(const qgstring &);
   qgstring* get();
   int isEmpty() const;
   int GetCount();
 private:
   int count;
   qgstring* UpperPtr; // pointer to the last added
                          // string in the stack
   qgstring* LowerPtr; // pointer to the first added
   
   qgstring* getNewqgs(const qgstring &);
};


int qgstack::GetCount()
{ return count;
};


/* constructor */
qgstack::qgstack()
{ UpperPtr=LowerPtr=0;
};

/* destructor */
qgstack::~qgstack()
{
   if (!isEmpty()) //stack not empty
     {
	qgstring* CurrentPtr=UpperPtr, *TempPtr;
	while(CurrentPtr!=0)
	  { TempPtr=CurrentPtr;
	     CurrentPtr=CurrentPtr->NextPtr;
	     delete TempPtr;
	  };       	
     };
   std::cout << "No more strings" << std::endl;
};

void qgstack::put(const qgstring & copy)
{
   count++;
   qgstring* NewPtr= getNewqgs(copy);
   if (isEmpty()) // stack empty
     UpperPtr=LowerPtr=NewPtr;
   else
     { //stack not empty
	NewPtr->NextPtr=UpperPtr;
	UpperPtr=NewPtr;
     };
};


qgstring* qgstack::get()
{
if (isEmpty()) return 0;
   else
     {	
   qgstring* TempPtr=UpperPtr;
   UpperPtr=TempPtr->NextPtr; // reducing the buffer
	count--;
//	std::cout << count << std::endl;
   return TempPtr; // and returning the element which we subtracted
     };   
};


int qgstack::isEmpty() const
{
   return UpperPtr == 0;
}

qgstring* qgstack::getNewqgs(const qgstring & copy)
{
   qgstring* ptr = new qgstring(copy);
   assert(ptr != 0);
   return ptr;
};
#endif