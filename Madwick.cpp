float qt[4] = { 0, };

void Madgwick( k4a_imu_sample_t* pImuSample, float gain )
{
	const float gx = pImuSample->gyro_sample.xyz.x * 0.5f;
	const float gy = pImuSample->gyro_sample.xyz.y * 0.5f;
	const float gz = pImuSample->gyro_sample.xyz.z * 0.5f;

	float qdp[4];
	// 角速度
	qdp[0] = -qt[1] * gx - qt[2] * gy - qt[3] * gz;
	qdp[1] =  qt[0] * gx - qt[3] * gy + qt[2] * gz;
	qdp[2] =  qt[3] * gx + qt[0] * gy - qt[1] * gz;
	qdp[3] = -qt[2] * gx + qt[1] * gy + qt[0] * gz;

	// 偶然のフィーバーを避ける
	if ( ! ((pImuSample->acc_sample.xyz.x == 0.0f) && (pImuSample->acc_sample.xyz.y == 0.0f) && (pImuSample->acc_sample.xyz.z == 0.0f)) )
	{
		// 加速度
		float fn = 1.0f / sqrtf( pImuSample->acc_sample.xyz.x * pImuSample->acc_sample.xyz.x + pImuSample->acc_sample.xyz.y * pImuSample->acc_sample.xyz.y + pImuSample->acc_sample.xyz.z * pImuSample->acc_sample.xyz.z );
		const float ax = pImuSample->acc_sample.xyz.x * fn;
		const float ay = pImuSample->acc_sample.xyz.y * fn;
		const float az = pImuSample->acc_sample.xyz.z * fn;

		float q2[4], q4[4], qn[4], qw[4], qx[4];
		for( int i = 0; i < 4; i++ )
		{
			q2[i] = qt[i] * 2.0f;
			q4[i] = qt[i] * 4.0f;
			qn[i] = qt[i] * 8.0f;
			qw[i] = qt[i] * qt[i];
			qx[i] = qw[i] * 4.0f;
		}

		qn[0] = q4[0] * qw[2] + q2[2] * ax + q4[0] * qw[1] - q2[1] * ay;
		qn[1] = q4[1] * qw[3] - q2[3] * ax + qt[1] * qx[0] - q2[0] * ay - q4[1] + qn[1] * qw[1] + qn[1] * qw[2] + q4[1] * az;
		qn[2] = qt[2] * qx[0] + q2[0] * ax + q4[2] * qw[3] - q2[3] * ay - q4[2] + qn[2] * q2[1] + qn[2] * q2[2] + q4[2] * az;
		qn[3] = qt[3] * qx[1] - q2[1] * ax + qt[3] * qx[2] - q2[2] * ay;
		fn = 1.0f / sqrtf( qn[0] * qn[0] + qn[1] * qn[1] + qn[2] * qn[2] + qn[3] * qn[3] ) * gain;
		qdp[0] = qdp[0] - qn[0] * fn;
		qdp[1] = qdp[1] - qn[1] * fn;
		qdp[2] = qdp[2] - qn[2] * fn;
		qdp[3] = qdp[3] - qn[3] * fn;
	}

	// IMU 周期ぶん積分してクォータニオンを正規化
	qt[0] += qdp[0] / 1600.0f;
	qt[1] += qdp[1] / 1600.0f;
	qt[2] += qdp[2] / 1600.0f;
	qt[3] += qdp[3] / 1600.0f;
	const float fn = 1.0f / sqrtf( qt[0] * qt[0] + qt[1] * qt[1] + qt[2] * qt[2] + qt[3] * qt[3] );
	qt[0] *= fn;
	qt[1] *= fn;
	qt[2] *= fn;
	qt[3] *= fn;
}
