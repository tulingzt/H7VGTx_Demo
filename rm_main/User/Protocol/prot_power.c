#include "prot_power.h"
#include "can_comm.h"
#include "prot_judge.h"
#include "string.h"

supercap_t supercap;
power_control_t power_control;

void power_init(void)
{
    memset(&supercap, 0, sizeof(supercap_t));
    memset(&power_control, 0, sizeof(power_control_t));
    power_control.judge_power_buffer = 60.0f;
    power_control.judge_max_power    = 55;
    power_control.min_buffer         = 30;
    power_control.limit_kp           = 0.4f;
    supercap.max_volage              = 23.6f;
    supercap.min_volage              = 13.5f;
}

void power_judge_update(void)
{
    power_control.judge_chassis_power = power_heat_data.chassis_power;
    power_control.judge_max_power     = robot_status.chassis_power_limit;
    power_control.judge_power_buffer  = power_heat_data.buffer_energy;
}

void power_output_data(void)
{
    uint8_t send_buff[8] = {0};
    send_buff[0] = supercap.mode;
    memcpy(send_buff + 1, &supercap.charge_current_ref, sizeof(supercap.charge_current_ref));
    can_std_transmit(CAN_CHANNEL_2, 0x021, send_buff);
}

void power_get_data(uint8_t *data)
{
    uint16_t supercap_voltage = 0;
    uint16_t source_power = 0;
    uint16_t chassis_power = 0;
    
    memcpy(&supercap_voltage, data, 2);
    memcpy(&source_power, data + 2, 2);
    memcpy(&chassis_power, data + 4, 2);
    //电容电压信息
    supercap.volage = supercap_voltage / 100.0f;
    supercap.volume_percent = (supercap.volage - supercap.min_volage) / (supercap.max_volage - supercap.min_volage) * 100.0f;
    //功率信息
    power_control.source_power = source_power / 100.0f;
    power_control.chassis_power = chassis_power / 100.0f;
}