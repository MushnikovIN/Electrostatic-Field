function validate_geometry(elementary_wires)
    %VALIDATE_GEOMETRY Проверка взаимного расположения элементарных проводов.
    %   VALIDATE_GEOMETRY(ELEMENTARY_WIRES) выполняет попарный расчет расстояний
    %   между всеми элементарными проводниками для обнаружения пересечений,
    %   а также контролирует отсутствие касания или пересечения поверхности земли.
    %
    %   При обнаружении геометрических наложений функция генерирует ошибку
    %   и прерывает выполнение программы.
    %
    %   Входные аргументы:
    %       ELEMENTARY_WIRES - Массив структур элементарных проводов                    

    % Определение количества элементарных проводников
    num_wires = length(elementary_wires);
    
    % проверка на касание или пересечение земли
    for i = 1:num_wires
        % Расчет расстояния от нижней точки провода до заземляющей плоскости (y = 0)
        distance_to_ground = elementary_wires(i).y - elementary_wires(i).r0;
        
        if distance_to_ground < 0
            error('Геометрическая ошибка: Элементарный провод №%d пересекает или касается земли! (y = %.4f, r0 = %.4f, нижняя точка = %.4f).', ...
                i, elementary_wires(i).y, elementary_wires(i).r0, distance_to_ground);
        end
    end
    
    % проверка на взаимное перекрытие проводов
    % Попарное сравнение всех уникальных пар проводников
    for i = 1:num_wires
        for j = (i + 1):num_wires
            % Расчет расстояния между центрами двух окружностей
            dist_centers = sqrt((elementary_wires(i).x - elementary_wires(j).x)^2 + (elementary_wires(i).y  - elementary_wires(j).y)^2);
            
            % Вычисление минимально допустимого расстояния между центрами (суммы радиусов)
            sum_radii = elementary_wires(i).r0 + elementary_wires(j).r0;
            
            % Генерация ошибки при обнаружении пересечения окружностей проводников
            if dist_centers < sum_radii
                error(['Геометрическая ошибка: Обнаружено пересечение проводов!\n' ...
                       'Элементарный провод №%d (x=%.4f, y=%.4f, r0=%.4f) перекрывается с\n' ...
                       'элементарным проводом №%d (x=%.4f, y=%.4f, r0=%.4f).\n' ...
                       'Расстояние между центрами: %.4f (минимум должно быть: %.4f).'], ...
                       i, elementary_wires(i).x, elementary_wires(i).y, elementary_wires(i).r0, ...
                       j, elementary_wires(j).x, elementary_wires(j).y, elementary_wires(j).r0, ...
                       dist_centers, sum_radii);
            end
        end
    end
    
    % Вывод уведомления об успешном прохождении тестов в командное окно
    fprintf('Геометрическая валидация успешна: пересечений и касаний земли не обнаружено.\n');
end
