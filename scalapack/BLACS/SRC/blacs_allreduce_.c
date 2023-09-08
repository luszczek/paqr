#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_allreduce(Int ConTxt, Int *msg)
#else
F_VOID_FUNC blacs_allreduce_(Int *ConTxt, Int *msg)
#endif
{
   extern BLACSCONTEXT **BI_MyContxts;
   extern Int BI_MaxNCtxt;
   BLACSCONTEXT *ctxt;
   Int ierr;
/*
 * Make sure context handle is in range
 */
   if ( (Mpval(ConTxt) >= 0) && (Mpval(ConTxt) < BI_MaxNCtxt) )
   {
/*
 *    Make sure context is still defined
 */
      ctxt = BI_MyContxts[Mpval(ConTxt)];
      if (ctxt != NULL)
      {
          long long msg2 = *msg;
          long long msg3 = 0;
          MPI_Request request;
          MPI_Status status;
//printf("~ %d %d %d\n", *msg, ctxt->cscp.Iam, ctxt->rscp.Iam);
//         ierr=MPI_Allreduce(MPI_IN_PLACE, msg, 1, MPI_LONG_LONG_INT, MPI_MAX, ctxt->scp->comm);
         ierr=MPI_Allreduce(&msg2, &msg3, 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);//ctxt->scp->comm);
//         ierr=MPI_Iallreduce(&msg2, &msg3, 1, MPI_LONG_LONG_INT, MPI_MAX, ctxt->scp->comm, &request);
//printf("! %d\n", ierr);
//         ierr = MPI_Wait(&request, &status);
//printf("@ %d\n", ierr);
//printf("@ %d %d\n", *msg, ierr);
         *msg = msg3;
//printf("# %d %d\n", *msg, ierr);
      }
   }
}
