% Очистка рабочей области и командного окна
clear
clc

% Задание исходной конфигурации фаз и параметров проводов системы
Wires = struct('x', {-4, 4, 10}, ...
    'y',            {4, 8, 4}, ...
    'U',            {1000, 500, 900}, ...
    'is_splited',   {false, false, true}, ...
    'r0',           {0.01, 0.02, 0.049}, ...
    'N',            {0, 0, 50}, ...
    'd',            {0, 0, 0.1});

% Первичная верификация исходной структуры
validate_wires(Wires)

% Преобразование расщепленных фаз в массив элементарных проводников
Wires = splited_to_elementary(Wires);

% Проверка геометрии элементарных проводов на пересечения и касание земли
validate_geometry(Wires)

% Задание относительной диэлектрической проницаемости среды
eps = 1;

% Инициализация объекта класса Field
field = Field(Wires, eps);

% Отрисовка координатной сетки и конфигурации проводов
field.plot_wires();

% Вывод тестовой расчетной точки на график для проверки работоспособности методов
field.plot_E_phi_at_point([1 1], "both")

% Инициализация флага продолжения итерационного цикла расчетов
need_culc = true;

% Запуск основного цикла интерактивного расчета параметров поля в точках
while need_culc
    % Запрос режима ввода координат
    graph1_or_keybord2 = input("Задать точки на графике или с клавиатуры? Введите 1 или 2 (1 - на графике; 2 - с клавиатуры):\n");
    
    if graph1_or_keybord2 == 1 || graph1_or_keybord2 == 2
        
        % Запрос отображаемого параметра
        correct_param = false;
        while ~correct_param
            param_choice = input("Какие параметры считать? Введите число:\n 1 - Только напряженность (E)\n 2 - Только потенциал (phi)\n 3 - И то, и другое (both)\n");
            switch param_choice
                case 1, plot_var = 'E'; correct_param = true;
                case 2, plot_var = 'phi'; correct_param = true;
                case 3, plot_var = 'both'; correct_param = true;
                otherwise
                    fprintf("Некорректный выбор параметра. Пожалуйста, введите 1, 2 или 3.\n\n");
            end
        end
        
        % Запрос количества рассчитываемых точек
        N = input("Сколько точек вы хотели бы рассчитать?\n");
        
        % Цикл обработки точек в зависимости от выбранного режима
        if graph1_or_keybord2 == 1
            % Ввод кликом по графику
            for idx = 1:N
                fprintf('Ожидание клика для точки %d из %d...\n', idx, N);
                
                % Интерактивное получение координат клика мыши
                coords = field.get_point_from_fig(); 
                
                % Отрисовка маркера и подписи физических величин на графике
                field.plot_E_phi_at_point(coords, plot_var); 
                
                % Вывод базовых координат точки в командное окно
                fprintf('Точка %d (график): X = %.4f, Y = %.4f.', idx, coords(1), coords(2));

                % расчет и вывод параметров поля в консоль
                switch plot_var
                    case 'E'
                        fprintf('  Расчетное значение: E = %.4f В/м\n\n', field.E_calculation(coords(1), coords(2)));
                    case 'phi'
                        fprintf('  Расчетное значение: phi = %.4f В\n\n', field.Phi_calculation(coords(1), coords(2)));
                    case 'both'
                        fprintf('  Расчетные значения: E = %.4f В/м, phi = %.4f В\n\n', field.E_calculation(coords(1), coords(2)), field.Phi_calculation(coords(1), coords(2)));
                end
            end
            
        else
            % Ввод вручную с клавиатуры 
            for idx = 1:N
                fprintf('Ввод данных для точки %d из %d:\n', idx, N);
                
                % Ручной ввод координат пользователем
                x_in = input('  Введите координату X: ');
                y_in = input('  Введите координату Y: ');
                
                % Формирование вектора координат и его визуализация
                coords = [x_in, y_in];
                field.plot_E_phi_at_point(coords, plot_var);
                
                % Вывод базовых координат точки в командное окно
                fprintf('Точка %d (клавиатура): X = %.4f, Y = %.4f.', idx, coords(1), coords(2));
                
                % Матричный расчет и динамический вывод физических параметров в консоль
                switch plot_var
                    case 'E'
                        fprintf('  Расчетное значение: E = %.4f В/м\n\n', field.E_calculation(coords(1), coords(2)));
                    case 'phi'
                        fprintf('  Расчетное значение: phi = %.4f В\n\n', field.Phi_calculation(coords(1), coords(2)));
                    case 'both'
                        fprintf('  Расчетные значения: E = %.4f В/м, phi = %.4f В\n\n', field.E_calculation(coords(1), coords(2)), field.Phi_calculation(coords(1), coords(2)));
                end
            end
        end
        
        % Запрос на повторение или завершение цикла расчетов
        user_continue = input("Желаете продолжить расчеты для новых точек? (Y - Да / N - Нет):\n", 's');
        if ~(strcmpi(user_continue, 'Y') || strcmpi(user_continue, 'yes') || strcmpi(user_continue, 'да'))
            need_culc = false;
            fprintf("Расчет окончен. Программа завершена.\n");
        end
        
    else
        % Обработка некорректного выбора режима ввода
        fprintf("Введено некорректное значение режима ввода. Повторите попытку.\n\n");
    end
end
