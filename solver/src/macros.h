
#ifndef NDEBUG
#define DEBUG(command) command;
#else
#define DEBUG(command)
#endif

#define EMPTY_ID -1 // empty field ID
#define P_ID 1      // field ID for pressure data fields
#define RHS_ID 2    // field ID for rhs data fields
#define U_ID 3      // field ID for u data fields
#define F_ID 4      // field ID for f data fields
#define V_ID 5      // field ID for v data fields
#define G_ID 6      // field ID for g data fields
#define Q_ID 7      // field ID for q data fields
