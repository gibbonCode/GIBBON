function plan = buildfile
import matlab.buildtool.tasks.*

plan = buildplan(localfunctions);

fc = matlab.buildtool.io.FileCollection.fromPaths(["docs/*.m", "lib/*.m", "tests/*.m"]);
plan("check") = CodeIssuesTask(fc);
plan.DefaultTasks = "check";
end
