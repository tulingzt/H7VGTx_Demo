#include "debug_task.h"
#include "cmsis_os.h"
#include "data_log.h"
#include "stdint.h"

#include "gimbal_task.h"
#include "shoot_task.h"
#include "wlr.h"
#include "leg_vmc.h"
#include "wheel_leg_model.h"
#include "drv_dji_motor.h"
#include "prot_judge.h"
#include "prot_power.h"
#include "prot_imu.h"
#include "kalman_filter.h"
#include "func_generator.h"
#include "prot_dr16.h"
#include "prot_judge.h"
#include "us_time.h"

us_time_t test_time;
kalman_filter_t test;
uint8_t debug_wave = 9;

void log_scope_data_pkg(void)
{
    switch(debug_wave) {
        case 1: {//云台pid调试
//            log_scope_get_data(gimbal.yaw_spd.ref);
//            log_scope_get_data(gimbal.yaw_spd.fdb);
//            log_scope_get_data(gimbal.yaw_angle.ref);
//            log_scope_get_data(gimbal.yaw_angle.fdb);
//            log_scope_get_data(gimbal.yaw_output);
//            log_scope_get_data(yaw_motor.tx_current);
            
            log_scope_get_data(gimbal.pit_spd.ref);
            log_scope_get_data(gimbal.pit_spd.fdb);
            log_scope_get_data(gimbal.pit_angle.ref);
            log_scope_get_data(gimbal.pit_angle.fdb);
//            log_scope_get_data(gimbal.pit_output);
//            log_scope_get_data(pit_motor.tx_current);
            break;
        } case 2: {//拨盘pid调试
//            log_scope_get_data(shoot.trigger_spd.ref);
//            log_scope_get_data(shoot.trigger_spd.fdb);
//            log_scope_get_data(shoot.trigger_ecd.ref);
//            log_scope_get_data(shoot.trigger_ecd.fdb);
//            log_scope_get_data(shoot.trigger_output);
//            log_scope_get_data(trigger_motor.tx_current);
            break;
        } case 3: {//底盘yaw roll调试
//            log_scope_get_data(wlr.wz_set);
//            log_scope_get_data(wlr.wz_fdb);
//            log_scope_get_data(wlr.yaw_fdb + wlr.yaw_err);
//            log_scope_get_data(wlr.yaw_fdb);
            
//            log_scope_get_data(wlr.roll_set);
//            log_scope_get_data(wlr.roll_fdb);
            break;
        } case 4: {//底盘功率调试
            log_scope_get_data(power_heat_data.chassis_power);
            log_scope_get_data(power_heat_data.buffer_energy);
            log_scope_get_data(supercap.volage);
            break;
        } case 5: {//支撑力调试
            log_scope_get_data(wlr.side[0].Fn_fdb);
            log_scope_get_data(wlr.side[0].Fn_kal);
            log_scope_get_data(wlr.side[0].fly_cnt);
            log_scope_get_data(wlr.side[1].Fn_fdb);
            log_scope_get_data(wlr.side[1].Fn_kal);
            log_scope_get_data(wlr.side[1].fly_cnt);
            break;
        } case 6: {//向心力调试
//            log_scope_get_data(twm.v0_ref);
            log_scope_get_data(twm.v0_fdb);
//            log_scope_get_data(twm.w0_ref);
            log_scope_get_data(twm.w0_imu_fdb);
//            log_scope_get_data(twm.a_ref);
            log_scope_get_data(twm.a_fdb);
//            log_scope_get_data(twm.r_ref);
            log_scope_get_data(twm.r_fdb);
            log_scope_get_data(twm.gravity_compensate_angle);
        } case 7: {
            log_scope_get_data(vmc[0].L_fdb);
            log_scope_get_data(vmc[1].L_fdb);
            log_scope_get_data(wlr.high_set);
        } case 8: {
            log_scope_get_data(wlr.side[0].Fy);
            log_scope_get_data(wlr.side[0].T0);
            log_scope_get_data(wlr.side[1].Fy);
            log_scope_get_data(wlr.side[1].T0);
        } case 9: {
            log_scope_get_data(wlr.side[0].v_fdb);
            log_scope_get_data(wlr.side[0].v_kal);
//            log_scope_get_data(wlr.side[1].v_fdb);
//            log_scope_get_data(wlr.side[1].v_kal);
            log_scope_get_data(wlr.side[0].a_fdb);
            log_scope_get_data(wlr.side[0].a_kal);
//            log_scope_get_data(wlr.side[1].a_fdb);
//            log_scope_get_data(wlr.side[1].a_kal);
        }default:break;
    }
}

/* 串口上位机数据发送任务 */
void debug_task(void const* argument)
{
    uint32_t thread_wake_time = osKernelSysTick();
    for(;;)
    {
        log_scope_data_output();
        osDelayUntil(&thread_wake_time, 5);
    }
}
