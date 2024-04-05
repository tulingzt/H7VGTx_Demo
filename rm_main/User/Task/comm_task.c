#include "comm_task.h"
#include "cmsis_os.h"
#include "prot_vision.h"
#include "drv_dji_motor.h"
#include "drv_ht_motor.h"

uint32_t temp_cnt;
void comm_task(void const* argument)
{
    uint32_t thread_wake_time = osKernelSysTick();
    for(;;)
    {
        temp_cnt++;
        dji_motor_output_data();
        ht_motor_output_data();
        vision_output_data();
//        ht_motor_output_single_data(&joint_motor[0]);
//        ht_motor_output_single_data(&joint_motor[1]);
        osDelayUntil(&thread_wake_time, 2);
//        dji_motor_output_data();
//        ht_motor_output_single_data(&joint_motor[2]);
//        ht_motor_output_single_data(&joint_motor[3]);
//        osDelayUntil(&thread_wake_time, 2);
    }
}
