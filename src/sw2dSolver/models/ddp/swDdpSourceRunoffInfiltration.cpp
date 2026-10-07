#include <cmath>
#include <vector>

#include "./swDdp.h"
#include "../../../common/swIO.h"
#include "../../../common/swLog.h"
#include "../common/swSource.h"
#include "../../../common/param.h"


/** @todo: adapt code */
void swDdp::sourceInfiltration(void)
/**
 * @brief advances the solution in time for the infiltration step
 * @return nothing
 * @par Modified attributes:
 *   - water_depth
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{
}


/** @todo: adapt code */
void swDdp::sourceRunoffEverywhere(void)
/**
 * @brief advances the solution in time for the runoff generation step when rain falls onto the complete cell
 * @return nothing
 * @par Modified attributes:
 *   - water_depth
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{
}


/** @todo: adapt code */
void swDdp::sourceRunoffWateronly(void)
/**
 * @brief advances the solution in time for the runoff generation step when rain falls onto water alone
 * @return nothing
 * @par Modified attributes:
 *   - water_depth
 *   - unit_discharge_x
 *   - unit_discharge_y
 */
{
}