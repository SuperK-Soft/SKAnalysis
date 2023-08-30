# SkipEventFlags
A simple Tool to reject events based on the event flags. Specify a set of allowed flags with:

	allowedFlags 1,2,3

to only allow events where event flag bits 1, 2 or 3 are set and reject any events where none of these bits are set. Alternatively specify a set of rejected  flags with:

	skippedFlags 11,12,13

to skip events where flag bits 11, 12 or 13 are set. Flag bits are retrieved from skhead_.ifevsk.

Triggers may be specified by their bit number, or by their name:
+ Flags for SK-I-III:
- "ATM" = 0
- "TRG" = 1
- "SMP_REGISTER" = 2
- "SCALER" = 3
- "PEDESTAL_START" = 4
- "PEDESTAL_DATA(ATM)" = 5
- "PEDESTAL_HISTOGRAM" = 6
- "PEDESTAL_END" = 7
- "END_OF_RUN" = 8
- "PEDESTAL_ON" = 9
- "unknown_10" = 10
- "GPS_DATA" = 11
- "CAMAC_ADC" = 12
- "ANTI_DATA" = 13
- "INNER_SLOW_DATA" = 14
- "RUN_INFORMATION" = 15
- "ERROR_TKO-PS" = 16
- "ERROR_HV-PS" = 17
- "ERROR_TEMPERATURE" = 18
- "unknown_19" = 19
- "UNCOMPLETED_ATM_DATA" = 20
- "INVALID_ATM_DATA" = 21
- "unknown_22" = 22
- "unknown_23" = 23
- "ERROR_DATA" = 24
- "UNREASONABLE_DATA" = 25
- "LED_BURST_ON" = 26
- "unknown_27" = 27
- "INNER_DETECTOR_OFF" = 28
- "ANTI_DETECTOR_OFF" = 29
- "unknown_30" = 30
- "TRG_IS_AVAILABLE" = 31

* Flags for SK-IV+
- "QBEE_TQ" = 0
- "HARD_TRG" = 1
- "QBEE_STAT" = 2
- "DB_STAT_BLOCK" = 3
- "CORRUPTED_CHECKSUM" = 4
- "MISSING_SPACER" = 5
- "PED_HIST_BLOCK" = 6
- "unknown_7" = 7
- "unknown_8" = 8
- "PEDESTAL_ON" = 9
- "RAW_AMT_BLOCK" = 10
- "GPS_DATA" = 11
- "PEDESTAL_CHECK" = 12
- "SEND_BLOCK" = 13
- "INNER_SLOW_DATA" = 14
- "RUN_INFORMATION" = 15
- "PREV_T0_BLOCK" = 16
- "unknown_17" = 17
- "FE_TRL_BLOCK" = 18
- "SPACER_BLOCK" = 19
- "INCOMPLETE_TQ" = 20
- "CORRUPT_TQ_BLOCK" = 21
- "TRG_MISMATCH_TQ" = 22
- "QBEE_ERROR" = 23
- "SORT_BLOCK" = 24
- "CORRUPTED_BLOCK" = 25
- "LED_BURST_ON" = 26
- "EVNT_TRAILER" = 27
- "INNER_DETECTOR_OFF" = 28
- "ANTI_DETECTOR_OFF" = 29
- "T2K_GPS" = 30
- "EVNT_HDR_&_SOFTWARE_TRG" = 31

