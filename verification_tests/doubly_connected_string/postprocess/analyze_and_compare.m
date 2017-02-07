combined = load('combined.txt');
combined = combined';
fcom = fft(combined);
%plot(abs(fcom(1:100,:)));


absfcom = abs(fcom);
onePlot = absfcom(:,5);
cutPlot = onePlot(1:100);
maxPlot = max(cutPlot);

for i = 1:100
	if (cutPlot(i) == maxPlot)
		maxFreq = i;
	end
end

tolerance = 1.0e-15;

refFreq = load('frequency.ref');

fprintf('refFreq and maxFreq = %d, %d\n', refFreq, maxFreq);

if (abs(refFreq - maxFreq) > tolerance)

	error('frequencies different');

end
