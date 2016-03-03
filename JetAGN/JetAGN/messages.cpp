/// Messages
#define MAX_BUFFER_SIZE 255
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "messages.h"

void show_message(MessageType message, ProgramModule module)
{
    show_time();
    static clock_t startTime[255];

    switch(message)
    {
        case msgStart:
            startTime[module] = clock();
            printf("Starting ");
            show_module(module);
            break;

        case msgEnd:
            printf("Finished ");
            show_module(module);
            printf("(%gs) ", diff_seconds(startTime[module], clock()) );
            break;

        case msgError:
            printf("ERROR in: ");
            show_module(module);
            break;
        
//        default:
  //          show_error(Invalid_MessageType, Module_Message);
    //        break;
    }
    printf("\n");    
    return;
}



// Imprime la hora %H:%M:%S en pantalla.
void show_time(void)
{
	
    time_t rawTime;
  //  struct tm *timeInfo;
    char buffer[MAX_BUFFER_SIZE];
    
	struct tm timeInfo;


	time(&rawTime);
	localtime_s(&timeInfo, &rawTime);


 //   time(&rawTime);     
 //   timeInfo = localtime(&rawTime);
    
    strftime(buffer, MAX_BUFFER_SIZE, "%H:%M:%S ", &timeInfo);
    printf("%s", buffer);
    
    return;
}

void show_module(ProgramModule module)
{

    switch(module)
    {
        case Module_Main:
            printf("main ");
            break;
	
        case Module_electronInjection:
            printf("Injection of electrons");
            break;
      
        case Module_Message:
            printf("message ");
            break;

		case Module_electronDistribution:
            printf("Calculating Electron Distribution ");
            break;


		case Module_luminosities:
			printf("Calculating non-thermal SEDs");
			break;

 
  /*      default:
            printf("\n");
            show_error(Invalid_Module, Module_Message);
            break;*/ 
	}

//    return;
}

double diff_seconds(clock_t timeStart, clock_t timeEnd)
{
    return ((double )(timeEnd - timeStart))/CLOCKS_PER_SEC;
}

/* void show_error(ErrorType error, ProgramModule module)
{
    printf("\t ");

    switch(error)
    {
        case Invalid_Energy:
            show_message(msgError, module);
            printf("Invalid energy.");
            break;

        case Invalid_Particle:
            show_message(msgError, module);
            printf("Invalid particle type.");
            break;

        case Invalid_Pointer:
            show_message(msgError, module);
            printf("Invalid pointer.");
            break;

        case Invalid_MessageType:
            show_message(msgError, module);
            printf("Invalid message type.");
            break;

        case Invalid_Module:
            show_message(msgError, Module_Message);
            printf("Invalid Module.");
            break;

        case Invalid_Process:
            show_message(msgError, module);
            printf("Invalid process.");
            break;

        default:
            show_message(msgError, Module_Message);
            printf("Unknown error message.");
            break;
    }

    printf("\n");
    system("pause");
    exit(1);
}*/ 